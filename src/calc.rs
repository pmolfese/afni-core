// PUBLIC DOMAIN NOTICE
//
// This file is part of afni-core, which was written by employees of the United
// States Government (National Institutes of Health) as part of their official
// duties. It is a "United States Government Work" (17 U.S.C. 105) and is in the
// public domain; outside the US, rights are waived under CC0 1.0. See LICENSE.
//
// ---------------------------------------------------------------------------
// WHAT THIS FILE IS
//
// An evaluator for `3dcalc`-style expressions: the language of `3dcalc`,
// `1deval` and `ccalc`, as implemented by AFNI's `parser.f`. It turns text such
// as `step(a-3)*step(b-2)` into a small program once (`Expr::parse`) and then
// evaluates it any number of times (`Expr::eval`), for example once per voxel
// to build a mask.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * Pure and file-neutral: it knows nothing about datasets or voxels. A caller
//   binds the one-letter variables `a`..`z` to numbers for each evaluation
//   (a dataset sub-brick value, a voxel coordinate, another mask).
// * Written for viewers that threshold and combine overlays ("mask where A and
//   B"), and reusable by sumaru for surface masks.
//
// FIDELITY (read parser.f, then compared with `ccalc`/`1deval`)
//
//   Grammar. Spaces are removed and letters upper-cased first, so `2 3` is 23
//   and `Step(A)` is `STEP(A)`. Precedence, lowest to highest: `+ -`, then
//   `* /`, then `**` (also written `^`), which associates to the RIGHT
//   (`2^3^2` = 512). A unary minus applies to the whole power expression that
//   follows it (`-2^2` = -4, `2^-1` = 0.5); a unary plus is ignored. `[ ]` are
//   the same as `( )`. AFNI has NO relational, boolean or conditional
//   operators (no `>`, `&&`, `?:`, `%`); its masks are built with `step`,
//   `astep`, `within`, `equals`, `and`, `or`, `not` and `ifelse`. Those all
//   work here exactly as in 3dcalc, and this crate ADDS C-style operators as
//   sugar (see DEPARTURES): `< <= > >= == !=`, `&& || !`, and `c ? t : f`.
//
//   "Designed not to fail" (parser.f line 496): illegal operations become
//   legal ones instead of NaN: `x/0` is 0, `sqrt(x)` is `sqrt(|x|)`, `log(x)`
//   is `log(|x|)` (and 0 for 0), `asin`/`acos`/`atanh`/`acosh` leave an
//   out-of-domain argument unchanged, `exp` caps its argument at 87.5, and
//   `a**b` returns `a` unchanged unless `a > 0` or (`a != 0` and `b` is an
//   integer). These are reproduced exactly.
//
//   Comparisons follow Fortran/C: with a NaN argument `step(x)` is 1 (`x <= 0`
//   is false), as in AFNI. Callers that hold missing data should decide what a
//   missing value is before binding it (afniru binds NaN as 0).
//
// DEPARTURES (also listed in docs/DIFFERENCES_FROM_AFNI.md)
//
//   * ADDED operators, absent from AFNI (which rejects them). Precedence,
//     lowest to highest: `?:` (right-associative), `||`, `&&`, `== !=`,
//     `< <= > >=`, `+ -`, `* /`, unary `! - +`, `**`. Comparisons and boolean
//     operators give 1 or 0; `&&`, `||`, `!` and the `?` test treat any
//     nonzero value as true (like `and`, `or`, `not`, `ifelse`). Unlike
//     `step`, a comparison with NaN is false. `a<b` is `step(b-a)`, `a==b` is
//     `equals(a,b)`, `a&&b` is `and(a,b)`, `c?t:f` is `ifelse(c,t,f)`; the
//     conformance tests check these identities against AFNI itself. A single
//     `=`, `&` or `|` is an error rather than being guessed at.
//   * Only the core function set is implemented (see `FUNCTIONS`). The rest of
//     AFNI's list (random numbers, Bessel/Airy/gamma/erf, the 27 `fico_*`
//     conversions, `hrfbk*`, `rhddc2`, `acfwxm`, `gamp/gamq`, `isprime`,
//     `ztone`, `cdf2stat/stat2cdf`) is recognised by name and rejected with an
//     `Unsupported` error naming the function, never evaluated approximately.
//   * Variables are single letters, as in 3dcalc; any other identifier is an
//     error naming it (parser.f accepts longer symbols that 3dcalc then
//     rejects).
//   * All whitespace is ignored, not only the space character.
//   * `absextreme` returns the largest absolute value, as AFNI documents. AFNI's
//     own evaluator never reaches that code (an 8-character opcode is compared
//     with the 10-character name) and returns the argument count instead.
//   * `within` needs 3 arguments, `orstat` 2 or more, `pairmax`/`pairmin` 2 or
//     more (AFNI reads past the arguments for fewer).
// ---------------------------------------------------------------------------

//! `3dcalc`-style expressions: parse once, evaluate many times.
//!
//! ```
//! use afni_core::calc::Expr;
//!
//! // Voxels where `a` exceeds 3 and `b` exceeds 2, as 3dcalc would write it.
//! let both = Expr::parse("step(a-3)*step(b-2)")?;
//! assert_eq!(both.variables(), ['a', 'b']);
//! assert_eq!(both.eval(|v| if v == 'a' { 4.0 } else { 2.5 }), 1.0);
//! assert_eq!(both.eval(|v| if v == 'a' { 4.0 } else { 1.0 }), 0.0);
//! # Ok::<(), afni_core::Error>(())
//! ```

use crate::error::{Error, Result};

/// Degrees to radians, `PI / 180`.
const D2R: f64 = std::f64::consts::PI / 180.0;

/// One instruction of a compiled expression.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Op {
    /// Push a constant.
    Num(f64),
    /// Push the value of variable `a`..`z` (index 0..26).
    Var(u8),
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Neg,
    // Extensions (not in AFNI): comparisons and boolean logic, all 1.0 / 0.0.
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
    And,
    Or,
    Not,
    /// `c ? t : f`: pops three values, pushes `t` if `c` is not zero, else `f`.
    Select,
    /// Call a function with this many arguments (already on the stack).
    Call(Func, u8),
}

macro_rules! functions {
    ($( $variant:ident => $name:literal , $args:expr ;)*) => {
        /// The functions this evaluator implements.
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        #[allow(missing_docs)]
        enum Func { $($variant),* }

        /// The implemented functions with their argument counts (`None` means
        /// "one or more", AFNI's variable-argument functions). Names are
        /// upper case, as AFNI folds the expression to upper case.
        pub const FUNCTIONS: &[(&str, Option<usize>)] = &[ $(($name, $args)),* ];

        impl Func {
            fn from_name(name: &str) -> Option<Func> {
                match name { $($name => Some(Func::$variant),)* _ => None }
            }
        }
    };
}

functions! {
    Sin => "SIN", Some(1); Cos => "COS", Some(1); Tan => "TAN", Some(1);
    Asin => "ASIN", Some(1); Acos => "ACOS", Some(1); Atan => "ATAN", Some(1);
    Atan2 => "ATAN2", Some(2);
    Sinh => "SINH", Some(1); Cosh => "COSH", Some(1); Tanh => "TANH", Some(1);
    Asinh => "ASINH", Some(1); Acosh => "ACOSH", Some(1); Atanh => "ATANH", Some(1);
    Exp => "EXP", Some(1); Log => "LOG", Some(1); Log10 => "LOG10", Some(1);
    Abs => "ABS", Some(1); Int => "INT", Some(1); Sqrt => "SQRT", Some(1);
    Cbrt => "CBRT", Some(1);
    Max => "MAX", Some(2); Min => "MIN", Some(2); Mod => "MOD", Some(2);
    Sind => "SIND", Some(1); Cosd => "COSD", Some(1); Tand => "TAND", Some(1);
    Rect => "RECT", Some(1); Step => "STEP", Some(1); Bool => "BOOL", Some(1);
    Posval => "POSVAL", Some(1); Tent => "TENT", Some(1); Bell2 => "BELL2", Some(1);
    Notzero => "NOTZERO", Some(1); Iszero => "ISZERO", Some(1); Not => "NOT", Some(1);
    Ispositive => "ISPOSITIVE", Some(1); Isnegative => "ISNEGATIVE", Some(1);
    Equals => "EQUALS", Some(2); Astep => "ASTEP", Some(2);
    Ifelse => "IFELSE", Some(3);
    And => "AND", None; Or => "OR", None; Mofn => "MOFN", None;
    Within => "WITHIN", None; Amongst => "AMONGST", None;
    Median => "MEDIAN", None; Mean => "MEAN", None; Stdev => "STDEV", None;
    Sem => "SEM", None; Mad => "MAD", None; Orstat => "ORSTAT", None;
    Argmax => "ARGMAX", None; Argnum => "ARGNUM", None; Choose => "CHOOSE", None;
    Pairmax => "PAIRMAX", None; Pairmin => "PAIRMIN", None;
    Minabove => "MINABOVE", None; Maxbelow => "MAXBELOW", None;
    Extreme => "EXTREME", None; Absextreme => "ABSEXTREME", None;
    Lmode => "LMODE", None; Hmode => "HMODE", None;
}

/// Functions AFNI has that this evaluator rejects by name.
const UNSUPPORTED: &[&str] = &[
    "AI", "DAI", "I0", "I1", "J0", "J1", "K0", "K1", "Y0", "Y1", "BI", "DBI", "ERF", "ERFC",
    "GAMMA", "QG", "QGINV", "GRAN", "URAN", "IRAN", "ERAN", "LRAN", "PLEG", "ZTONE", "CDF2STAT",
    "STAT2CDF", "RHDDC2", "HRFBK4", "HRFBK5", "LOGCOSH", "ACFWXM", "GAMP", "GAMQ", "ISPRIME",
    "FICO_T2P", "FICO_P2T", "FICO_T2Z", "FITT_T2P", "FITT_P2T", "FITT_T2Z", "FIFT_T2P", "FIFT_P2T",
    "FIFT_T2Z", "FIZT_T2P", "FIZT_P2T", "FIZT_T2Z", "FICT_T2P", "FICT_P2T", "FICT_T2Z", "FIBT_T2P",
    "FIBT_P2T", "FIBT_T2Z", "FIBN_T2P", "FIBN_P2T", "FIBN_T2Z", "FIGT_T2P", "FIGT_P2T", "FIGT_T2Z",
    "FIPT_T2P", "FIPT_P2T", "FIPT_T2Z",
];

fn bad(reason: impl Into<String>) -> Error {
    Error::InvalidParameter {
        name: "expression".into(),
        reason: reason.into(),
    }
}

// ---------------------------------------------------------------------------
// Tokens
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Number(f64),
    /// A name that is not a function (a variable, or an error later).
    Name(String),
    Func(Func),
    Plus,
    Minus,
    Star,
    Slash,
    Power,
    Open,
    Close,
    Comma,
    // Extensions (not in AFNI).
    Less,
    LessEq,
    Greater,
    GreaterEq,
    EqEq,
    NotEq,
    AndAnd,
    OrOr,
    Bang,
    Question,
    Colon,
    End,
}

/// Split the (space-stripped, upper-cased) text into tokens, as `GET_TOKEN`.
fn tokenize(text: &str) -> Result<Vec<Token>> {
    let chars: Vec<char> = text.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        match c {
            '+' => tokens.push(Token::Plus),
            '-' => tokens.push(Token::Minus),
            '/' => tokens.push(Token::Slash),
            '*' if chars.get(i + 1) == Some(&'*') => {
                tokens.push(Token::Power);
                i += 1;
            }
            '*' => tokens.push(Token::Star),
            '^' => tokens.push(Token::Power),
            '(' | '[' => tokens.push(Token::Open),
            ')' | ']' => tokens.push(Token::Close),
            ',' => tokens.push(Token::Comma),
            '<' | '>' | '!' | '=' | '&' | '|' => {
                let two = chars.get(i + 1) == Some(&'=');
                let token = match (c, two) {
                    ('<', true) => Token::LessEq,
                    ('<', false) => Token::Less,
                    ('>', true) => Token::GreaterEq,
                    ('>', false) => Token::Greater,
                    ('!', true) => Token::NotEq,
                    ('!', false) => Token::Bang,
                    ('=', true) => Token::EqEq,
                    ('=', false) => {
                        return Err(bad(
                            "a single '=' is not an operator: write '==' to compare",
                        ));
                    }
                    ('&', _) | ('|', _) => {
                        if chars.get(i + 1) != Some(&c) {
                            return Err(bad(format!(
                                "write '{c}{c}' for logical {}",
                                if c == '&' { "and" } else { "or" }
                            )));
                        }
                        if c == '&' {
                            Token::AndAnd
                        } else {
                            Token::OrOr
                        }
                    }
                    _ => unreachable!(),
                };
                if matches!(token, Token::AndAnd | Token::OrOr) || two {
                    i += 1;
                }
                tokens.push(token);
            }
            '?' => tokens.push(Token::Question),
            ':' => tokens.push(Token::Colon),
            c if c.is_ascii_alphabetic() => {
                let start = i;
                while i + 1 < chars.len()
                    && (chars[i + 1].is_ascii_alphanumeric()
                        || chars[i + 1] == '_'
                        || chars[i + 1] == '$')
                {
                    i += 1;
                }
                let name: String = chars[start..=i].iter().collect();
                if let Some(f) = Func::from_name(&name) {
                    tokens.push(Token::Func(f));
                } else if UNSUPPORTED.contains(&name.as_str()) {
                    return Err(Error::Unsupported(format!(
                        "the function {} is not implemented (supported: {})",
                        name.to_lowercase(),
                        FUNCTIONS
                            .iter()
                            .map(|(n, _)| n.to_lowercase())
                            .collect::<Vec<_>>()
                            .join(", ")
                    )));
                } else if name == "PI" {
                    tokens.push(Token::Number(std::f64::consts::PI));
                } else {
                    tokens.push(Token::Name(name));
                }
            }
            c if c.is_ascii_digit() || c == '.' => {
                let start = i;
                let digit = |k: usize| chars.get(k).is_some_and(char::is_ascii_digit);
                while digit(i + 1) {
                    i += 1;
                }
                // A decimal point, unless the number began with one.
                if c != '.' && chars.get(i + 1) == Some(&'.') {
                    i += 1;
                    while digit(i + 1) {
                        i += 1;
                    }
                }
                // An exponent counts only if a digit follows it.
                if matches!(chars.get(i + 1), Some('E') | Some('D')) {
                    let mut k = i + 2;
                    if matches!(chars.get(k), Some('+') | Some('-')) {
                        k += 1;
                    }
                    if digit(k) {
                        i = k;
                        while digit(i + 1) {
                            i += 1;
                        }
                    }
                }
                let literal: String = chars[start..=i]
                    .iter()
                    .collect::<String>()
                    .replace('D', "E");
                let value = literal
                    .parse::<f64>()
                    .map_err(|_| bad(format!("cannot read the number '{literal}'")))?;
                tokens.push(Token::Number(value));
            }
            other => return Err(bad(format!("cannot interpret the symbol '{other}'"))),
        }
        i += 1;
    }
    tokens.push(Token::End);
    Ok(tokens)
}

// ---------------------------------------------------------------------------
// Parser (recursive descent for AFNI's LL(1) grammar)
// ---------------------------------------------------------------------------

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
    code: Vec<Op>,
    used: [bool; 26],
}

impl Parser {
    fn peek(&self) -> &Token {
        &self.tokens[self.pos]
    }

    fn next(&mut self) -> Token {
        let t = self.tokens[self.pos].clone();
        if self.pos + 1 < self.tokens.len() {
            self.pos += 1;
        }
        t
    }

    /// Extension: `c ? t : f`, right-associative, lowest precedence.
    fn expr(&mut self) -> Result<()> {
        self.or()?;
        if *self.peek() == Token::Question {
            self.next();
            self.expr()?;
            if self.next() != Token::Colon {
                return Err(bad("expected a ':' to finish the '?' choice"));
            }
            self.expr()?;
            self.code.push(Op::Select);
        }
        Ok(())
    }

    /// Extension: `||`, left to right.
    fn or(&mut self) -> Result<()> {
        self.and()?;
        while *self.peek() == Token::OrOr {
            self.next();
            self.and()?;
            self.code.push(Op::Or);
        }
        Ok(())
    }

    /// Extension: `&&`, left to right.
    fn and(&mut self) -> Result<()> {
        self.equality()?;
        while *self.peek() == Token::AndAnd {
            self.next();
            self.equality()?;
            self.code.push(Op::And);
        }
        Ok(())
    }

    /// Extension: `==` and `!=`.
    fn equality(&mut self) -> Result<()> {
        self.relation()?;
        loop {
            let op = match self.peek() {
                Token::EqEq => Op::Eq,
                Token::NotEq => Op::Ne,
                _ => return Ok(()),
            };
            self.next();
            self.relation()?;
            self.code.push(op);
        }
    }

    /// Extension: `<`, `<=`, `>`, `>=`.
    fn relation(&mut self) -> Result<()> {
        self.sum()?;
        loop {
            let op = match self.peek() {
                Token::Less => Op::Lt,
                Token::LessEq => Op::Le,
                Token::Greater => Op::Gt,
                Token::GreaterEq => Op::Ge,
                _ => return Ok(()),
            };
            self.next();
            self.sum()?;
            self.code.push(op);
        }
    }

    /// `E4`: sum of terms, left to right.
    fn sum(&mut self) -> Result<()> {
        self.term()?;
        loop {
            let op = match self.peek() {
                Token::Plus => Op::Add,
                Token::Minus => Op::Sub,
                _ => return Ok(()),
            };
            self.next();
            self.term()?;
            self.code.push(op);
        }
    }

    /// `E6`: product of factors, left to right.
    fn term(&mut self) -> Result<()> {
        self.power_expr()?;
        loop {
            let op = match self.peek() {
                Token::Star => Op::Mul,
                Token::Slash => Op::Div,
                _ => return Ok(()),
            };
            self.next();
            self.power_expr()?;
            self.code.push(op);
        }
    }

    /// `E9 E8`: an atom, then a right-associative power chain.
    fn power_expr(&mut self) -> Result<()> {
        self.atom()?;
        if *self.peek() == Token::Power {
            self.next();
            self.power_expr()?;
            self.code.push(Op::Pow);
        }
        Ok(())
    }

    /// `E9`: a number, variable, call, parenthesised sum, or a sign.
    fn atom(&mut self) -> Result<()> {
        match self.next() {
            Token::Number(v) => self.code.push(Op::Num(v)),
            Token::Name(name) => {
                let mut chars = name.chars();
                match (chars.next(), chars.next()) {
                    (Some(c), None) if c.is_ascii_uppercase() => {
                        let index = c as u8 - b'A';
                        self.used[index as usize] = true;
                        self.code.push(Op::Var(index));
                    }
                    _ => {
                        return Err(bad(format!(
                            "unknown symbol '{}': variables are single letters a to z",
                            name.to_lowercase()
                        )));
                    }
                }
            }
            Token::Func(f) => self.call(f)?,
            Token::Open => {
                self.expr()?;
                self.expect_close()?;
            }
            // A unary plus is dropped.
            Token::Plus => self.atom()?,
            // A unary minus applies to the atom and its power chain.
            Token::Minus => {
                self.power_expr()?;
                self.code.push(Op::Neg);
            }
            // A unary `!` binds like a unary minus.
            Token::Bang => {
                self.power_expr()?;
                self.code.push(Op::Not);
            }
            Token::Star
            | Token::Slash
            | Token::Power
            | Token::Comma
            | Token::Less
            | Token::LessEq
            | Token::Greater
            | Token::GreaterEq
            | Token::EqEq
            | Token::NotEq
            | Token::AndAnd
            | Token::OrOr
            | Token::Question
            | Token::Colon => {
                return Err(bad("expected a value before an operator"));
            }
            Token::Close => return Err(bad("expected an expression before ')'")),
            Token::End => return Err(bad("unexpected end of the expression")),
        }
        Ok(())
    }

    fn expect_close(&mut self) -> Result<()> {
        match self.next() {
            Token::Close => Ok(()),
            Token::End => Err(bad("expected a \")\" before the end of the expression")),
            _ => Err(bad("expected an operator or a \")\"")),
        }
    }

    fn call(&mut self, f: Func) -> Result<()> {
        let name = FUNCTIONS
            .iter()
            .find(|(n, _)| Func::from_name(n) == Some(f))
            .map_or("?", |(n, _)| *n);
        if self.next() != Token::Open {
            return Err(bad(format!(
                "expected a \"(\" after {}",
                name.to_lowercase()
            )));
        }
        let mut count = 1usize;
        self.expr()?;
        while *self.peek() == Token::Comma {
            self.next();
            self.expr()?;
            count += 1;
        }
        self.expect_close()?;
        let wanted = FUNCTIONS
            .iter()
            .find(|(n, _)| *n == name)
            .and_then(|(_, a)| *a);
        if let Some(n) = wanted {
            if n != count {
                return Err(bad(format!(
                    "{} takes {n} argument{}, not {count}",
                    name.to_lowercase(),
                    if n == 1 { "" } else { "s" }
                )));
            }
        }
        let minimum = match f {
            Func::Within => 3,
            Func::Orstat | Func::Pairmax | Func::Pairmin => 2,
            _ => 1,
        };
        if count < minimum {
            return Err(bad(format!(
                "{} needs at least {minimum} arguments",
                name.to_lowercase()
            )));
        }
        self.code.push(Op::Call(f, count as u8));
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// The compiled expression
// ---------------------------------------------------------------------------

/// A parsed expression, ready to evaluate.
#[derive(Debug, Clone, PartialEq)]
pub struct Expr {
    text: String,
    code: Vec<Op>,
    used: [bool; 26],
}

impl Expr {
    /// Parse `text`. Errors say what was wrong (and, for functions AFNI has
    /// but this crate does not, which function).
    pub fn parse(text: &str) -> Result<Expr> {
        let cleaned: String = text
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .collect();
        if cleaned.is_empty() {
            return Err(bad("the expression is empty"));
        }
        if cleaned.len() > 9999 {
            return Err(bad("the expression is too long"));
        }
        let mut p = Parser {
            tokens: tokenize(&cleaned)?,
            pos: 0,
            code: Vec::new(),
            used: [false; 26],
        };
        p.expr()?;
        match p.peek() {
            Token::End => {}
            Token::Close => return Err(bad("unbalanced \")\"")),
            _ => return Err(bad("expected an operator")),
        }
        Ok(Expr {
            text: text.trim().to_string(),
            code: p.code,
            used: p.used,
        })
    }

    /// The expression as given.
    pub fn text(&self) -> &str {
        &self.text
    }

    /// The variables the expression uses, as lower-case letters in
    /// alphabetical order.
    pub fn variables(&self) -> Vec<char> {
        (0..26u8)
            .filter(|i| self.used[*i as usize])
            .map(|i| (b'a' + i) as char)
            .collect()
    }

    /// Evaluate with `value(letter)` giving each variable (lower-case letter)
    /// its number. Variables the expression does not use are never asked for.
    pub fn eval(&self, value: impl Fn(char) -> f64) -> f64 {
        let mut vars = [0.0; 26];
        for (i, v) in vars.iter_mut().enumerate() {
            if self.used[i] {
                *v = value((b'a' + i as u8) as char);
            }
        }
        self.eval_vars(&vars, &mut Vec::new())
    }

    /// Evaluate with all 26 variables given (`vars[0]` is `a`), reusing
    /// `stack` between calls so that evaluating millions of voxels does not
    /// allocate.
    pub fn eval_vars(&self, vars: &[f64; 26], stack: &mut Vec<f64>) -> f64 {
        stack.clear();
        for op in &self.code {
            match *op {
                Op::Num(v) => stack.push(v),
                Op::Var(i) => stack.push(vars[i as usize]),
                Op::Add => binary(stack, |a, b| a + b),
                Op::Sub => binary(stack, |a, b| a - b),
                Op::Mul => binary(stack, |a, b| a * b),
                Op::Div => binary(stack, |a, b| if b != 0.0 { a / b } else { 0.0 }),
                Op::Pow => binary(stack, power),
                Op::Neg => {
                    if let Some(top) = stack.last_mut() {
                        *top = -*top;
                    }
                }
                Op::Lt => binary(stack, |a, b| flag(a < b)),
                Op::Le => binary(stack, |a, b| flag(a <= b)),
                Op::Gt => binary(stack, |a, b| flag(a > b)),
                Op::Ge => binary(stack, |a, b| flag(a >= b)),
                Op::Eq => binary(stack, |a, b| flag(a == b)),
                Op::Ne => binary(stack, |a, b| flag(a != b)),
                Op::And => binary(stack, |a, b| flag(a != 0.0 && b != 0.0)),
                Op::Or => binary(stack, |a, b| flag(a != 0.0 || b != 0.0)),
                Op::Not => {
                    if let Some(top) = stack.last_mut() {
                        *top = flag(*top == 0.0);
                    }
                }
                Op::Select => {
                    let f = stack.pop().unwrap_or(0.0);
                    let t = stack.pop().unwrap_or(0.0);
                    if let Some(c) = stack.last_mut() {
                        *c = if *c != 0.0 { t } else { f };
                    }
                }
                Op::Call(f, n) => {
                    let at = stack.len() - n as usize;
                    let result = call(f, &mut stack[at..]);
                    stack.truncate(at);
                    stack.push(result);
                }
            }
        }
        stack.pop().unwrap_or(0.0)
    }
}

fn binary(stack: &mut Vec<f64>, f: impl Fn(f64, f64) -> f64) {
    let b = stack.pop().unwrap_or(0.0);
    if let Some(a) = stack.last_mut() {
        *a = f(*a, b);
    }
}

/// AFNI's `**`: the left operand is returned unchanged unless the power is
/// well defined for it.
fn power(a: f64, b: f64) -> f64 {
    if a > 0.0 || (a != 0.0 && b == b.trunc()) {
        a.powf(b)
    } else {
        a
    }
}

fn step(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else {
        1.0
    }
}

fn flag(b: bool) -> f64 {
    if b {
        1.0
    } else {
        0.0
    }
}

fn boolean(x: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        1.0
    }
}

/// Sort ascending (NaN last), like the bubble sort of parser.f for ordinary data.
fn sorted(x: &[f64]) -> Vec<f64> {
    let mut v = x.to_vec();
    v.sort_by(f64::total_cmp);
    v
}

fn median(x: &[f64]) -> f64 {
    match x.len() {
        0 => 0.0,
        1 => x[0],
        2 => 0.5 * (x[0] + x[1]),
        n => {
            let s = sorted(x);
            if n % 2 == 0 {
                0.5 * (s[n / 2 - 1] + s[n / 2])
            } else {
                s[n / 2]
            }
        }
    }
}

fn mean(x: &[f64]) -> f64 {
    match x.len() {
        1 => x[0],
        2 => 0.5 * (x[0] + x[1]),
        n => x.iter().sum::<f64>() / n as f64,
    }
}

fn stdev(x: &[f64]) -> f64 {
    let n = x.len();
    if n <= 1 {
        return 0.0;
    }
    let bar = x.iter().sum::<f64>() / n as f64;
    (x.iter().map(|v| (v - bar).powi(2)).sum::<f64>() / (n as f64 - 1.0)).sqrt()
}

/// The most common value; ties go to the lowest (`low`) or highest value.
fn mode(x: &[f64], low: bool) -> f64 {
    if x.len() == 1 {
        return x[0];
    }
    let s = sorted(x);
    let (mut value, mut run, mut best_count, mut best) = (s[0], 1usize, 0usize, s[0]);
    for &v in &s[1..] {
        if v != value {
            if (low && run > best_count) || (!low && run >= best_count) {
                best = value;
                best_count = run;
            }
            value = v;
            run = 1;
        } else {
            run += 1;
        }
    }
    if (low && run > best_count) || (!low && run >= best_count) {
        best = value;
    }
    best
}

/// Evaluate one function on its arguments (`a` may be reordered).
fn call(f: Func, a: &mut [f64]) -> f64 {
    let x = a[0];
    let n = a.len();
    match f {
        Func::Sin => x.sin(),
        Func::Cos => x.cos(),
        Func::Tan => x.tan(),
        Func::Sind => (D2R * x).sin(),
        Func::Cosd => (D2R * x).cos(),
        Func::Tand => (D2R * x).tan(),
        Func::Asin if x.abs() <= 1.0 => x.asin(),
        Func::Acos if x.abs() <= 1.0 => x.acos(),
        Func::Asin | Func::Acos => x,
        Func::Atan => x.atan(),
        Func::Atan2 if x != 0.0 || a[1] != 0.0 => x.atan2(a[1]),
        Func::Atan2 => x,
        Func::Sinh if x.abs() < 87.5 => x.sinh(),
        Func::Cosh if x.abs() < 87.5 => x.cosh(),
        Func::Sinh | Func::Cosh => x,
        Func::Tanh => x.tanh(),
        Func::Asinh => {
            let ax = x.abs();
            let y = if ax <= 10.0 {
                ax + (ax * ax + 1.0).sqrt()
            } else {
                ax * (1.0 + (1.0 + (1.0 / ax).powi(2)).sqrt())
            };
            if x < 0.0 {
                -y.ln()
            } else {
                y.ln()
            }
        }
        Func::Acosh if x >= 1.0 => {
            let y = if x <= 10.0 {
                x + (x * x - 1.0).sqrt()
            } else {
                x * (1.0 + (1.0 - (1.0 / x).powi(2)).sqrt())
            };
            y.ln()
        }
        Func::Acosh => x,
        Func::Atanh if x.abs() < 1.0 => 0.5 * ((1.0 + x) / (1.0 - x)).ln(),
        Func::Atanh => x,
        Func::Exp => x.min(87.5).exp(),
        Func::Log if x != 0.0 => x.abs().ln(),
        Func::Log10 if x != 0.0 => x.abs().log10(),
        Func::Log | Func::Log10 => x,
        Func::Abs => x.abs(),
        Func::Int => x.trunc(),
        Func::Sqrt => x.abs().sqrt(),
        Func::Cbrt => x.cbrt(),
        Func::Max => x.max(a[1]),
        Func::Min => x.min(a[1]),
        Func::Mod => {
            if a[1] != 0.0 {
                x - a[1] * (x / a[1]).trunc()
            } else {
                0.0
            }
        }
        Func::Rect => {
            if x.abs() <= 0.5 {
                1.0
            } else {
                0.0
            }
        }
        Func::Step | Func::Ispositive => step(x),
        Func::Isnegative => step(-x),
        Func::Bool | Func::Notzero => boolean(x),
        Func::Iszero | Func::Not => 1.0 - boolean(x),
        Func::Posval => {
            if x <= 0.0 {
                0.0
            } else {
                x
            }
        }
        Func::Tent => {
            if x.abs() >= 1.0 {
                0.0
            } else {
                1.0 - x.abs()
            }
        }
        Func::Bell2 => {
            let ax = x.abs();
            if ax <= 0.5 {
                1.0 - 1.333_333_333_333_333_3 * ax * ax
            } else if ax <= 1.5 {
                0.666_666_666_666_667 * (1.5 - ax).powi(2)
            } else {
                0.0
            }
        }
        Func::Equals => 1.0 - boolean(x - a[1]),
        Func::Astep => {
            if x.abs() > a[1] {
                1.0
            } else {
                0.0
            }
        }
        Func::Ifelse => {
            if x != 0.0 {
                a[1]
            } else {
                a[2]
            }
        }
        Func::And => boolean(if a.iter().all(|v| *v != 0.0) {
            1.0
        } else {
            0.0
        }),
        Func::Or => {
            if a.iter().any(|v| *v != 0.0) {
                1.0
            } else {
                0.0
            }
        }
        Func::Mofn => {
            let m = x as i64;
            let c = a[1..].iter().filter(|v| **v != 0.0).count() as i64;
            if c >= m {
                1.0
            } else {
                0.0
            }
        }
        Func::Within => {
            if a[0] < a[1] || a[0] > a[2] {
                0.0
            } else {
                1.0
            }
        }
        Func::Amongst => {
            if a[1..].contains(&x) {
                1.0
            } else {
                0.0
            }
        }
        Func::Median => median(a),
        Func::Mean => mean(a),
        Func::Stdev => stdev(a),
        Func::Sem => stdev(a) / (n as f64 + 0.000001).sqrt(),
        Func::Mad => {
            if n == 1 {
                0.0
            } else if n == 2 {
                0.5 * (a[0] - a[1]).abs()
            } else {
                let m = median(a);
                let deviations: Vec<f64> = a.iter().map(|v| (v - m).abs()).collect();
                median(&deviations)
            }
        }
        Func::Orstat => {
            let rest = &a[1..];
            if rest.len() <= 1 {
                rest[0]
            } else {
                let i = (x as i64).clamp(1, rest.len() as i64) as usize;
                sorted(rest)[i - 1]
            }
        }
        Func::Argmax => {
            let (mut top, mut at, mut zeros) = (a[0], 1usize, usize::from(a[0] == 0.0));
            for (i, v) in a.iter().enumerate().skip(1) {
                if *v > top {
                    at = i + 1;
                    top = *v;
                }
                if *v == 0.0 {
                    zeros += 1;
                }
            }
            if zeros == n {
                0.0
            } else {
                at as f64
            }
        }
        Func::Argnum => a.iter().filter(|v| **v != 0.0).count() as f64,
        Func::Choose => {
            let m = x as i64;
            let count = n as i64 - 1;
            if m < 1 || count < m {
                0.0
            } else {
                a[m as usize]
            }
        }
        Func::Pairmax | Func::Pairmin => {
            if n <= 2 {
                return a[1];
            }
            let m = n / 2;
            let (mut best, mut pair) = (a[0], a[m]);
            for i in 1..m {
                let better = if f == Func::Pairmax {
                    a[i] > best
                } else {
                    a[i] < best
                };
                if better {
                    best = a[i];
                    pair = a[m + i];
                }
            }
            pair
        }
        Func::Minabove => {
            if n == 1 {
                return x;
            }
            let b = a[1..]
                .iter()
                .copied()
                .filter(|v| *v > x && *v < 1.0e38)
                .fold(1.0e38, f64::min);
            if b == 1.0e38 {
                x
            } else {
                b
            }
        }
        Func::Maxbelow => {
            if n == 1 {
                return x;
            }
            let b = a[1..]
                .iter()
                .copied()
                .filter(|v| *v < x && *v > -1.0e38)
                .fold(-1.0e38, f64::max);
            if b == -1.0e38 {
                x
            } else {
                b
            }
        }
        Func::Extreme => {
            if n == 1 {
                return x;
            }
            // AFNI quirk, kept: each |value| is compared with the previous
            // winner's *signed* value, so `extreme(1,-5,3)` is 3, not -5.
            let mut best = 0.0_f64;
            for v in a.iter() {
                if v.abs() > best {
                    best = *v;
                }
            }
            if best == 0.0 {
                x
            } else {
                best
            }
        }
        Func::Absextreme => {
            if n == 1 {
                return x;
            }
            let best = a
                .iter()
                .fold(0.0_f64, |b, v| if v.abs() > b { v.abs() } else { b });
            if best == 0.0 {
                x
            } else {
                best
            }
        }
        Func::Lmode => mode(a, true),
        Func::Hmode => mode(a, false),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(text: &str) -> f64 {
        Expr::parse(text).unwrap().eval(|_| 0.0)
    }

    fn with(text: &str, vars: &[(char, f64)]) -> f64 {
        Expr::parse(text)
            .unwrap()
            .eval(|c| vars.iter().find(|(v, _)| *v == c).map_or(0.0, |(_, x)| *x))
    }

    fn err(text: &str) -> Error {
        Expr::parse(text).unwrap_err()
    }

    #[test]
    fn precedence_associativity_and_unary_minus_match_afni() {
        // Values from `ccalc -eval`.
        assert_eq!(ev("2+3*4"), 14.0);
        assert_eq!(ev("(2+3)*4"), 20.0);
        assert_eq!(ev("10/4"), 2.5);
        assert_eq!(ev("5-3-1"), 1.0);
        assert_eq!(ev("8/4/2"), 1.0);
        assert_eq!(ev("2^3^2"), 512.0); // right associative
        assert_eq!(ev("2**3**2"), 512.0);
        assert_eq!(ev("-2^2"), -4.0); // unary minus is below the power
        assert_eq!(ev("-2**2"), -4.0);
        assert_eq!(ev("2^-1"), 0.5);
        assert_eq!(ev("2*-3"), -6.0);
        assert_eq!(ev("+5"), 5.0);
        assert_eq!(ev("[2+3]*4"), 20.0); // brackets are parentheses
    }

    #[test]
    fn spaces_vanish_and_case_is_ignored() {
        assert_eq!(ev("2 3"), 23.0); // AFNI joins the digits
        assert_eq!(ev("Step ( 1 )"), 1.0);
        assert_eq!(ev("STEP(1)"), 1.0);
        assert!((ev("pi") - std::f64::consts::PI).abs() < 1e-15);
        assert_eq!(with("A+B", &[('a', 1.0), ('b', 2.0)]), 3.0);
    }

    #[test]
    fn numbers_are_read_like_afni() {
        assert_eq!(ev("1e3"), 1000.0);
        assert_eq!(ev("1.5e-2"), 0.015);
        assert_eq!(ev(".5"), 0.5);
        assert_eq!(ev("3."), 3.0);
        assert_eq!(ev("1d2"), 100.0); // Fortran D exponent
    }

    #[test]
    fn illegal_operations_become_legal_ones() {
        assert_eq!(ev("1/0"), 0.0);
        assert_eq!(ev("0/0"), 0.0);
        assert_eq!(ev("sqrt(-4)"), 2.0);
        assert_eq!(ev("log(0)"), 0.0);
        assert_eq!(ev("log(-1)"), 0.0); // log|x|
        assert_eq!(ev("acos(2)"), 2.0); // out of domain: unchanged
        assert_eq!(ev("asin(-3)"), -3.0);
        assert_eq!(ev("atanh(2)"), 2.0);
        assert_eq!(ev("acosh(0.5)"), 0.5);
        assert_eq!(ev("mod(7,0)"), 0.0);
        assert!((ev("exp(1000)") - 87.5_f64.exp()).abs() / 87.5_f64.exp() < 1e-12);
        assert_eq!(ev("(0-2)**0.5"), -2.0); // non-integer power of a negative: left operand
        assert_eq!(ev("0**2"), 0.0);
    }

    #[test]
    fn mask_functions_match_3dcalc() {
        assert_eq!(ev("step(0)"), 0.0);
        assert_eq!(ev("step(-1)"), 0.0);
        assert_eq!(ev("step(0.001)"), 1.0);
        assert_eq!(ev("astep(3,2)"), 1.0);
        assert_eq!(ev("astep(-3,2)"), 1.0);
        assert_eq!(ev("astep(2,2)"), 0.0);
        assert_eq!(ev("within(2,1,3)"), 1.0);
        assert_eq!(ev("within(3,1,3)"), 1.0);
        assert_eq!(ev("within(3.1,1,3)"), 0.0);
        assert_eq!(ev("rect(0.5)"), 1.0);
        assert_eq!(ev("rect(0.6)"), 0.0);
        assert_eq!(ev("bool(0)"), 0.0);
        assert_eq!(ev("bool(0.1)"), 1.0);
        assert_eq!(ev("not(0)"), 1.0);
        assert_eq!(ev("not(2)"), 0.0);
        assert_eq!(ev("equals(2,2)"), 1.0);
        assert_eq!(ev("equals(2,3)"), 0.0);
        assert_eq!(ev("ifelse(1,5,6)"), 5.0);
        assert_eq!(ev("ifelse(0,5,6)"), 6.0);
        assert_eq!(ev("and(1,1,0)"), 0.0);
        assert_eq!(ev("and(1,2,3)"), 1.0);
        assert_eq!(ev("or(0,0,0)"), 0.0);
        assert_eq!(ev("or(0,0,3)"), 1.0);
        assert_eq!(ev("mofn(2,1,0,1)"), 1.0);
        assert_eq!(ev("mofn(3,1,0,1)"), 0.0);
        assert_eq!(ev("posval(-2)"), 0.0);
        assert_eq!(ev("posval(2)"), 2.0);
        assert_eq!(ev("ispositive(0)"), 0.0);
        assert_eq!(ev("isnegative(-1)"), 1.0);
    }

    #[test]
    fn a_conjunction_of_two_thresholds_is_a_product_of_steps() {
        let both = Expr::parse("step(a-3)*step(b-2)").unwrap();
        assert_eq!(both.variables(), ['a', 'b']);
        let run = |a, b| both.eval(|c| if c == 'a' { a } else { b });
        assert_eq!(run(4.0, 3.0), 1.0);
        assert_eq!(run(3.0, 3.0), 0.0); // step(0) = 0: strictly greater
        assert_eq!(run(4.0, 2.0), 0.0);
    }

    #[test]
    fn math_functions_match_afni() {
        assert_eq!(ev("int(-2.7)"), -2.0);
        assert_eq!(ev("int(2.7)"), 2.0);
        assert_eq!(ev("mod(7,3)"), 1.0);
        assert_eq!(ev("mod(-7,3)"), -1.0);
        assert_eq!(ev("abs(-3)"), 3.0);
        assert_eq!(ev("min(3,2)"), 2.0);
        assert_eq!(ev("max(3,2)"), 3.0);
        assert!((ev("sind(30)") - 0.5).abs() < 1e-12);
        assert!((ev("cosd(60)") - 0.5).abs() < 1e-12);
        assert_eq!(ev("cbrt(-8)"), -2.0);
        assert!((ev("atan2(1,1)") - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
    }

    #[test]
    fn statistics_of_arguments_match_afni() {
        assert_eq!(ev("median(3,1,2)"), 2.0);
        assert_eq!(ev("median(4,1,3,2)"), 2.5);
        assert_eq!(ev("mean(1,2,6)"), 3.0);
        assert_eq!(ev("stdev(1,2,3)"), 1.0);
        assert_eq!(ev("argmax(1,3,2)"), 2.0);
        assert_eq!(ev("argmax(0,0)"), 0.0);
        assert_eq!(ev("argnum(1,0,2)"), 2.0);
        assert_eq!(ev("choose(2,10,20,30)"), 20.0);
        assert_eq!(ev("choose(4,10,20,30)"), 0.0);
        assert_eq!(ev("amongst(2,1,2,3)"), 1.0);
        assert_eq!(ev("amongst(5,1,2,3)"), 0.0);
        assert_eq!(ev("orstat(1,5,3,9)"), 3.0);
        assert_eq!(ev("pairmin(3,2,7,5,-1,-2,-3,-4)"), -2.0); // 3dcalc -help
        assert_eq!(ev("pairmax(1,5,2,10,50,20)"), 50.0);
        assert_eq!(ev("lmode(1,2,2,3,3)"), 2.0);
        assert_eq!(ev("hmode(1,2,2,3,3)"), 3.0);
        assert_eq!(ev("minabove(2,1,3,5)"), 3.0);
        assert_eq!(ev("maxbelow(4,1,3,5)"), 3.0);
        // AFNI compares |v| with the previous winner's signed value (a quirk).
        assert_eq!(ev("extreme(1,-5,3)"), 3.0);
        assert_eq!(ev("extreme(1,3,-5)"), -5.0);
        assert_eq!(ev("absextreme(1,-5,3)"), 5.0);
    }

    #[test]
    fn variables_are_the_single_letters_used() {
        let e = Expr::parse("a*b + c/a").unwrap();
        assert_eq!(e.variables(), ['a', 'b', 'c']);
        assert_eq!(Expr::parse("3+2").unwrap().variables(), Vec::<char>::new());
        // Only used variables are asked for.
        let asked = std::cell::RefCell::new(Vec::new());
        Expr::parse("b+d").unwrap().eval(|c| {
            asked.borrow_mut().push(c);
            1.0
        });
        assert_eq!(*asked.borrow(), ['b', 'd']);
    }

    fn run(text: &str, a: f64, b: f64) -> f64 {
        Expr::parse(text)
            .unwrap()
            .eval(|c| if c == 'a' { a } else { b })
    }

    #[test]
    fn comparisons_give_one_or_zero() {
        for (t, a, b, want) in [
            ("a<b", 1.0, 2.0, 1.0),
            ("a<b", 2.0, 2.0, 0.0),
            ("a<=b", 2.0, 2.0, 1.0),
            ("a>b", 3.0, 2.0, 1.0),
            ("a>b", 2.0, 2.0, 0.0),
            ("a>=b", 2.0, 2.0, 1.0),
            ("a==b", 2.0, 2.0, 1.0),
            ("a==b", 2.0, 3.0, 0.0),
            ("a!=b", 2.0, 3.0, 1.0),
            ("a!=b", 2.0, 2.0, 0.0),
            ("a>-3", -2.0, 0.0, 1.0),
            ("a > 3", 3.5, 0.0, 1.0),
        ] {
            assert_eq!(run(t, a, b), want, "{t} with a={a} b={b}");
        }
        // NaN compares false, never true.
        assert_eq!(run("a<b", f64::NAN, 1.0), 0.0);
        assert_eq!(run("a>=b", f64::NAN, 1.0), 0.0);
        assert_eq!(run("a!=b", f64::NAN, 1.0), 1.0);
    }

    #[test]
    fn boolean_operators_treat_nonzero_as_true() {
        assert_eq!(run("a&&b", 5.0, -2.0), 1.0);
        assert_eq!(run("a&&b", 5.0, 0.0), 0.0);
        assert_eq!(run("a||b", 0.0, 0.0), 0.0);
        assert_eq!(run("a||b", 0.0, 0.1), 1.0);
        assert_eq!(run("!a", 0.0, 0.0), 1.0);
        assert_eq!(run("!a", 7.0, 0.0), 0.0);
        assert_eq!(run("!!a", 7.0, 0.0), 1.0);
    }

    #[test]
    fn operator_precedence_runs_from_ternary_down_to_power() {
        // && binds tighter than ||.
        assert_eq!(run("1||0&&0", 0.0, 0.0), 1.0);
        // Relations bind tighter than equality: (1<2)==(2<3).
        assert_eq!(run("1<2==2<3", 0.0, 0.0), 1.0);
        // Arithmetic binds tighter than relations.
        assert_eq!(run("a+1>b*2", 3.0, 2.0), 0.0);
        assert_eq!(run("a+1>=b*2", 3.0, 2.0), 1.0);
        // ! is unary: it applies before == but after **.
        assert_eq!(run("!a==b", 0.0, 1.0), 1.0);
        assert_eq!(run("!2^2", 0.0, 0.0), 0.0);
        // Chained comparisons are C-style, left to right: (1<2)<3.
        assert_eq!(run("3>2>1", 0.0, 0.0), 0.0);
        // A comparison can be a factor of a product.
        assert_eq!(run("(a>1)*10+b", 2.0, 5.0), 15.0);
    }

    #[test]
    fn the_conditional_picks_a_branch_and_nests_to_the_right() {
        assert_eq!(run("a>0 ? 10 : 20", 1.0, 0.0), 10.0);
        assert_eq!(run("a>0 ? 10 : 20", -1.0, 0.0), 20.0);
        // a ? b : c ? d : e is a ? b : (c ? d : e).
        let t = "a<0 ? -1 : a==0 ? 0 : 1";
        assert_eq!(run(t, -5.0, 0.0), -1.0);
        assert_eq!(run(t, 0.0, 0.0), 0.0);
        assert_eq!(run(t, 5.0, 0.0), 1.0);
        // Inside a function argument and inside parentheses.
        assert_eq!(run("abs(a<0 ? a : -a)", 3.0, 0.0), 3.0);
        assert_eq!(run("1+(a ? 2 : 3)", 0.0, 0.0), 4.0);
        // Variables in every branch are reported.
        let e = Expr::parse("a ? b : c").unwrap();
        assert_eq!(e.variables(), vec!['a', 'b', 'c']);
    }

    #[test]
    fn operators_agree_with_the_function_forms() {
        let f = [-2.0, -0.5, 0.0, 0.5, 2.0];
        for a in f {
            for b in f {
                for (op, func) in [
                    ("a<b", "step(b-a)"),
                    ("a>b", "step(a-b)"),
                    ("a<=b", "1-step(a-b)"),
                    ("a>=b", "1-step(b-a)"),
                    ("a==b", "equals(a,b)"),
                    ("a!=b", "1-equals(a,b)"),
                    ("a&&b", "and(a,b)"),
                    ("a||b", "or(a,b)"),
                    ("!a", "not(a)"),
                    ("a?b:7", "ifelse(a,b,7)"),
                ] {
                    assert_eq!(run(op, a, b), run(func, a, b), "{op} vs {func} at {a},{b}");
                }
            }
        }
    }

    #[test]
    fn errors_name_the_problem() {
        let msg = |t: &str| err(t).to_string();
        assert!(msg("").contains("empty"));
        assert!(msg("3 %").contains("cannot interpret"), "{}", msg("3 %"));
        assert!(msg("a>").contains("end"));
        assert!(msg("a=3").contains("'=='"));
        assert!(msg("a&b").contains("'&&'"));
        assert!(msg("a|b").contains("'||'"));
        assert!(msg("a?1").contains("':'"));
        assert!(msg("<3").contains("before an operator"));
        assert!(msg("1:2").contains("expected an operator"));
        assert!(msg("2 +").contains("end"));
        assert!(msg("(2+3").contains("\")\""));
        assert!(msg("2+3)").contains("unbalanced"));
        assert!(msg("step()").contains("expression"));
        assert!(msg("step(1,2)").contains("takes 1 argument"));
        assert!(msg("min(1)").contains("takes 2 arguments"));
        assert!(msg("step+3").contains("\"(\""));
        assert!(msg("step 3").contains("unknown symbol 'step3'")); // spaces vanish
        assert!(msg("a b").contains("unknown symbol 'ab'"), "{}", msg("a b"));
        assert!(msg("foo").contains("unknown symbol 'foo'"));
        assert!(msg("2a").contains("expected an operator"));
        assert!(msg("*2").contains("before an operator"));
        assert!(msg("within(1,2)").contains("at least 3"));
    }

    #[test]
    fn functions_afni_has_but_we_lack_are_named_in_the_error() {
        for name in [
            "gran(1,2)",
            "erf(0.5)",
            "j0(1)",
            "fitt_t2p(2,10)",
            "isprime(7)",
            "hrfbk5(1,2)",
        ] {
            match err(name) {
                Error::Unsupported(m) => {
                    let f = name.split('(').next().unwrap();
                    assert!(m.contains(f), "{m}");
                }
                other => panic!("{name}: {other:?}"),
            }
        }
    }

    #[test]
    fn nan_arguments_follow_the_fortran_comparisons() {
        // As in AFNI: `x <= 0` is false for NaN, so step(NaN) is 1.
        assert_eq!(with("step(a)", &[('a', f64::NAN)]), 1.0);
        assert!(with("a+1", &[('a', f64::NAN)]).is_nan());
        assert_eq!(with("bool(a)", &[('a', f64::NAN)]), 1.0);
    }

    #[test]
    fn evaluation_reuses_its_stack() {
        let e = Expr::parse("max(a,b)*2").unwrap();
        let mut stack = Vec::new();
        let mut vars = [0.0; 26];
        for k in 0..100 {
            vars[0] = k as f64;
            vars[1] = 50.0;
            assert_eq!(e.eval_vars(&vars, &mut stack), 2.0 * f64::from(k.max(50)));
        }
        assert!(stack.capacity() < 16);
    }

    #[test]
    fn the_function_table_lists_every_implemented_function_once() {
        let mut names: Vec<_> = FUNCTIONS.iter().map(|(n, _)| *n).collect();
        let before = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), before);
        assert!(names.contains(&"STEP") && names.contains(&"IFELSE"));
        // Nothing is both implemented and listed as unsupported.
        assert!(UNSUPPORTED.iter().all(|u| !names.contains(u)));
    }
}
