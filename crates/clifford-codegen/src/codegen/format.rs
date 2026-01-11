//! Code formatting utilities.
//!
//! This module provides utilities for formatting generated Rust code
//! using rustfmt.

use proc_macro2::TokenStream;
use std::io::Write;
use std::process::{Command, Stdio};

/// Formats a token stream using rustfmt.
///
/// Returns the formatted code as a string. If rustfmt is not available
/// or fails, returns the unformatted token stream as a string.
///
/// # Example
///
/// ```
/// use quote::quote;
/// use clifford_codegen::codegen::format_tokens;
///
/// let tokens = quote! {
///     pub struct Foo { x: i32 }
/// };
///
/// let formatted = format_tokens(&tokens);
/// // formatted will be properly indented and spaced
/// ```
pub fn format_tokens(tokens: &TokenStream) -> String {
    let unformatted = tokens.to_string();
    format_code(&unformatted).unwrap_or(unformatted)
}

/// Formats Rust code using rustfmt.
///
/// Returns `None` if rustfmt is not available or fails.
fn format_code(code: &str) -> Option<String> {
    let mut child = Command::new("rustfmt")
        .arg("--edition=2024")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .ok()?;

    // Write code to stdin
    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(code.as_bytes()).ok()?;
    }

    // Read formatted output
    let output = child.wait_with_output().ok()?;

    if output.status.success() {
        String::from_utf8(output.stdout).ok()
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quote::quote;

    #[test]
    fn formats_simple_struct() {
        let tokens = quote! {
            pub struct Foo{x:i32,y:i32}
        };

        let formatted = format_tokens(&tokens);

        // Check that it compiles (basic structure is preserved)
        assert!(formatted.contains("struct Foo"));
        assert!(formatted.contains("x"));
        assert!(formatted.contains("y"));
    }

    #[test]
    fn handles_complex_code() {
        let tokens = quote! {
            impl<T: Float> Vector<T> {
                pub fn new(x: T, y: T) -> Self {
                    Self { x, y }
                }

                pub fn x(&self) -> T {
                    self.x
                }
            }
        };

        let formatted = format_tokens(&tokens);
        assert!(formatted.contains("impl<T: Float> Vector<T>"));
        assert!(formatted.contains("pub fn new"));
    }
}
