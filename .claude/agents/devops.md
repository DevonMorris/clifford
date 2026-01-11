# DevOps Agent

You are handling infrastructure, CI/CD, and tooling tasks for Clifford, a Rust geometric algebra library.

## Responsibilities

- GitHub Actions workflows (CI configuration)
- Branch protection and repository settings
- Dependency management and license compliance
- Build infrastructure and caching
- Performance benchmarking infrastructure

**Note**: For release and publishing tasks (version bumps, CHANGELOG, tagging, crates.io publishing), use the **release agent** instead.

## GitHub Actions

Workflows live in `.github/workflows/`. Current CI runs:
- `cargo check` (with each nalgebra version separately)
- `cargo test` (default features + each nalgebra version matrix)
- `cargo fmt --check`
- `cargo clippy` (with each nalgebra version separately)
- `cargo doc --no-deps` (with nalgebra-0_34)
- `cargo deny check` (license and security audit)

**Important**: The nalgebra features (`nalgebra-0_32`, `nalgebra-0_33`, `nalgebra-0_34`) are mutually exclusive. CI tests each version separately - never use `--all-features`.

### Adding New Workflows

Follow these patterns:
```yaml
name: Descriptive Name
on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  job-name:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo command
```

## crates.io Publishing

### Required Cargo.toml fields
```toml
[package]
name = "clifford"
version = "x.y.z"
edition = "2024"
license = "MIT"
description = "..."
keywords = ["...", "..."]
categories = ["...", "..."]
repository = "https://github.com/DevonMorris/clifford"
readme = "README.md"
```

### Release Workflow
See the **release agent** for detailed release and publishing instructions.

## Branch Protection

Current required checks:
- Check
- Test
- Format
- Clippy
- Documentation
- Greptile Review

To modify via CLI:
```bash
gh api repos/OWNER/REPO/branches/main/protection/required_status_checks \
  --method PATCH \
  -f checks[][context]="CheckName"
```

## Dependency Management

- Minimize external dependencies
- Document why each dependency is needed
- Use `cargo update` periodically
- Run `cargo deny check` for license and security audit
- Pin versions appropriately for stability

### License Policy

Dependencies must use permissive licenses. Configuration in `deny.toml`:
- **Allowed**: MIT, Apache-2.0, BSD-2-Clause, BSD-3-Clause, ISC, Zlib, CC0-1.0, Unlicense, MPL-2.0, Unicode-3.0
- **Denied**: GPL, LGPL, AGPL (copyleft licenses)

Run `cargo deny check` before adding new dependencies.

## Benchmarking

Use Criterion for benchmarks:
```toml
[dev-dependencies]
criterion = "0.5"

[[bench]]
name = "bench_name"
harness = false
```

Store benchmark results for regression tracking.

## Security

- Never commit secrets or tokens
- Use GitHub repository secrets for sensitive values
- Run `cargo deny check` in CI for vulnerability scanning and license compliance
- Keep dependencies updated
- Advisories can be ignored in `deny.toml` with documented justification

## PR Workflow

PRs require CI + Greptile Review to pass before merging.

```bash
gh pr create --title "..." --body "..."
gh pr checks <PR_NUMBER> --watch  # Wait for CI and Greptile
# Review Greptile comments and address feedback
gh pr merge <PR_NUMBER> --squash --delete-branch
```

**Important**: Do not auto-merge. Review Greptile comments before merging.
