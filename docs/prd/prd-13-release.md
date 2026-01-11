# PRD-13: Release and Publishing Process

**Status**: Draft
**Goal**: Establish a robust process for versioning, changelog management, and publishing to crates.io

## Motivation

As clifford matures toward production use, we need a clear release process that:

1. Communicates breaking changes clearly to users
2. Maintains a useful changelog for upgrade decisions
3. Automates publishing to reduce human error
4. Preserves backwards compatibility where possible
5. Provides predictable upgrade paths

## Versioning Strategy

### Semantic Versioning

We follow [Semantic Versioning 2.0.0](https://semver.org/):

- **MAJOR** (X.0.0): Breaking API changes
- **MINOR** (0.X.0): New features, backwards compatible
- **PATCH** (0.0.X): Bug fixes, backwards compatible

### Pre-1.0 Versioning

While at `0.x.y`, we treat minor versions as potentially breaking:

| Version | Meaning |
|---------|---------|
| `0.1.0` | Initial release, API unstable |
| `0.1.x` | Bug fixes only, API stable within 0.1 |
| `0.2.0` | May contain breaking changes from 0.1 |

**Recommendation**: Users should pin to `0.1` (not `0`) in Cargo.toml:

```toml
# Good: allows 0.1.x patches
clifford = "0.1"

# Risky: allows 0.2.0 which may break
clifford = "0"
```

### Post-1.0 Versioning

After 1.0.0, we commit to strict semver:

- Breaking changes require major version bump
- Deprecations allowed in minor versions, removals in major
- MSRV changes are breaking (major version bump)

### What Constitutes a Breaking Change

**Breaking** (requires major version bump):

- Removing or renaming public types, functions, methods, or modules
- Changing function signatures (parameters, return types)
- Changing struct fields (public structs)
- Changing trait bounds in ways that break existing impls
- Changing behavior in ways that violate documented contracts
- Removing trait implementations
- Increasing MSRV (Minimum Supported Rust Version)
- Changing default feature flags in ways that affect compilation

**Not Breaking** (minor version bump acceptable):

- Adding new public items (types, functions, modules)
- Adding new trait implementations (with care for coherence)
- Adding new optional features
- Deprecating items (but not removing)
- Bug fixes that match documented behavior
- Performance improvements
- Documentation changes

**Edge Cases**:

- Adding required trait methods: Breaking (use default impl instead)
- Adding fields to `#[non_exhaustive]` structs: Not breaking
- Changing private implementation details: Not breaking

## CHANGELOG Management

### Format

We use [Keep a Changelog](https://keepachangelog.com/) format:

```markdown
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- New `Motor::slerp()` method for smooth interpolation (#42)

### Changed
- Improved performance of geometric product by 15% (#45)

### Deprecated
- `Rotor::from_euler()` - use `Rotor::from_euler_angles()` instead (#43)

### Removed
- Removed deprecated `Vector::mag()` method (#41)

### Fixed
- Fixed normalization edge case for zero vectors (#44)

### Security
- Updated `foo` dependency to address CVE-XXXX-YYYY (#46)

## [0.2.0] - 2026-02-15

### Added
- Projective Geometric Algebra (PGA) support (#30)
- Rerun visualization integration (#39)

### Changed
- **BREAKING**: Renamed `dim3::Vec3` to `dim3::Vector` (#35)

## [0.1.0] - 2026-01-01

### Added
- Initial release with Euclidean 2D/3D support
- Generic multivector implementation
- nalgebra integration
```

### Changelog Rules

1. **Every PR updates the changelog** under `[Unreleased]`
2. **Use present tense**: "Add feature" not "Added feature"
3. **Reference PR/issue numbers**: `(#42)` at end of line
4. **Mark breaking changes**: Prefix with `**BREAKING**:`
5. **Group related changes**: One bullet can cover multiple related items
6. **Be user-focused**: Describe impact, not implementation details

### Release Process Updates Changelog

When releasing, the `[Unreleased]` section becomes the new version:

```markdown
## [Unreleased]

## [0.2.0] - 2026-02-15
<!-- Former [Unreleased] content moves here -->
```

## Release Workflow

### Pre-Release Checklist

Before creating a release:

1. **All CI passes on main**
2. **CHANGELOG.md updated**
   - `[Unreleased]` section has all changes
   - No placeholder text
3. **Version in Cargo.toml matches intended release**
4. **Documentation builds without warnings**
   ```bash
   cargo doc --no-deps --all-features
   ```
5. **All deprecation warnings addressed or documented**
6. **Breaking changes clearly documented in CHANGELOG**
7. **README reflects current API** (examples work)
8. **MSRV documented and tested**

### Release Steps

1. **Create release branch**
   ```bash
   git checkout main
   git pull origin main
   git checkout -b release/v0.2.0
   ```

2. **Update version**
   ```bash
   # Edit Cargo.toml: version = "0.2.0"
   cargo check  # Updates Cargo.lock
   ```

3. **Finalize CHANGELOG**
   - Move `[Unreleased]` content to new version section
   - Add release date
   - Add empty `[Unreleased]` section

4. **Create PR for release**
   ```bash
   git add Cargo.toml Cargo.lock CHANGELOG.md
   git commit -m "chore: release v0.2.0"
   git push origin release/v0.2.0
   gh pr create --title "Release v0.2.0" --body "Release preparation for v0.2.0"
   ```

5. **Merge PR after CI passes**

6. **Create and push tag**
   ```bash
   git checkout main
   git pull origin main
   git tag -a v0.2.0 -m "Release v0.2.0"
   git push origin v0.2.0
   ```

7. **Verify publish workflow succeeds**
   - Check GitHub Actions for publish job
   - Verify crate appears on crates.io

8. **Create GitHub Release**
   ```bash
   gh release create v0.2.0 --title "v0.2.0" --notes-file RELEASE_NOTES.md
   ```
   Or use the CHANGELOG section for that version as release notes.

### Automated Publishing

The existing `publish.yml` workflow handles publishing:

```yaml
on:
  push:
    tags:
      - 'v*.*.*'
```

**Enhancements needed**:

1. **Verify CHANGELOG has entry for version**
2. **Create GitHub Release automatically**
3. **Post to social media / announce** (optional)

## Backwards Compatibility

### Deprecation Policy

1. **Deprecate before removing**: Items must be deprecated for at least one minor version before removal
2. **Provide migration path**: Deprecation message must explain what to use instead
3. **Document in CHANGELOG**: All deprecations listed under `### Deprecated`

```rust
#[deprecated(
    since = "0.2.0",
    note = "Use `Rotor::from_euler_angles()` instead"
)]
pub fn from_euler(roll: T, pitch: T, yaw: T) -> Self {
    Self::from_euler_angles(roll, pitch, yaw)
}
```

### Feature Stability

| Feature | Stability | Policy |
|---------|-----------|--------|
| Core types (Vector, Rotor, etc.) | Stable | No breaking changes without major bump |
| Generic Multivector | Stable | No breaking changes without major bump |
| nalgebra integration | Stable per version | Each `nalgebra-0_XX` feature is stable |
| rerun integration | Unstable | May change in minor versions (0.x) |
| Visualization helpers | Unstable | May change in minor versions |

### MSRV Policy

- **Current MSRV**: Rust edition 2024 (1.85.0+)
- **MSRV bumps are breaking changes** (major version after 1.0)
- **Document MSRV** in Cargo.toml and README
- **Test MSRV in CI**

```toml
# Cargo.toml
[package]
rust-version = "1.85.0"
```

```yaml
# ci.yml addition
msrv:
  name: MSRV
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@1.85.0
    - run: cargo check
```

## Pre-Release Versions

For testing before official release:

### Alpha/Beta Releases

```toml
version = "0.2.0-alpha.1"
version = "0.2.0-beta.1"
version = "0.2.0-rc.1"
```

- **alpha**: Early testing, API may change significantly
- **beta**: Feature complete, API mostly stable
- **rc**: Release candidate, only bug fixes expected

### Publishing Pre-Releases

```bash
git tag -a v0.2.0-alpha.1 -m "v0.2.0-alpha.1"
git push origin v0.2.0-alpha.1
```

Pre-releases publish to crates.io but aren't selected by default:

```toml
# Gets 0.1.x, not 0.2.0-alpha.1
clifford = "0.1"

# Explicitly request pre-release
clifford = "0.2.0-alpha.1"
```

## CI Enhancements

### Release Validation Job

Add to `ci.yml`:

```yaml
release-ready:
  name: Release Ready
  runs-on: ubuntu-latest
  if: github.ref == 'refs/heads/main'
  steps:
    - uses: actions/checkout@v4

    - name: Check CHANGELOG has Unreleased section
      run: |
        if ! grep -q "## \[Unreleased\]" CHANGELOG.md; then
          echo "CHANGELOG.md missing [Unreleased] section"
          exit 1
        fi

    - name: Check version not already published
      run: |
        VERSION=$(cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')
        if cargo search clifford --limit 1 | grep -q "clifford = \"$VERSION\""; then
          echo "Version $VERSION already published - bump version for next release"
          exit 1
        fi
```

### Enhanced Publish Workflow

Update `publish.yml`:

```yaml
name: Publish to crates.io

on:
  push:
    tags:
      - 'v*.*.*'

env:
  CARGO_TERM_COLOR: always

jobs:
  publish:
    name: Publish
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: dtolnay/rust-toolchain@stable

      - uses: Swatinem/rust-cache@v2

      - name: Verify version matches tag
        run: |
          CARGO_VERSION=$(cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')
          TAG_VERSION=${GITHUB_REF#refs/tags/v}
          if [ "$CARGO_VERSION" != "$TAG_VERSION" ]; then
            echo "Version mismatch: Cargo.toml has $CARGO_VERSION, tag is $TAG_VERSION"
            exit 1
          fi

      - name: Verify CHANGELOG has release entry
        run: |
          TAG_VERSION=${GITHUB_REF#refs/tags/v}
          if ! grep -q "## \[$TAG_VERSION\]" CHANGELOG.md; then
            echo "CHANGELOG.md missing entry for version $TAG_VERSION"
            exit 1
          fi

      - name: Run full CI checks
        run: |
          cargo check
          cargo fmt --all -- --check
          cargo clippy -- -D warnings
          cargo test
          cargo doc --no-deps

      - name: Dry run publish
        run: cargo publish --dry-run

      - name: Publish to crates.io
        run: cargo publish
        env:
          CARGO_REGISTRY_TOKEN: ${{ secrets.CARGO_REGISTRY_TOKEN }}

      - name: Create GitHub Release
        uses: softprops/action-gh-release@v1
        with:
          body_path: RELEASE_NOTES.md
          generate_release_notes: true
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
```

## Files to Create/Modify

### New Files
- `CHANGELOG.md` - Changelog following Keep a Changelog format
- `docs/prd/prd-13-release.md` - This PRD

### Modified Files
- `Cargo.toml` - Add `rust-version` field for MSRV
- `.github/workflows/publish.yml` - Enhanced publish workflow
- `.github/workflows/ci.yml` - Add MSRV and release-ready jobs
- `README.md` - Document MSRV and versioning policy
- `CLAUDE.md` - Document release process

## Verification Checklist

- [ ] CHANGELOG.md exists with proper format
- [ ] Cargo.toml has `rust-version` field
- [ ] CI tests MSRV
- [ ] Publish workflow verifies CHANGELOG
- [ ] GitHub Release created automatically
- [ ] README documents versioning policy
- [ ] Deprecation examples work correctly

## Summary

This PRD establishes:

1. **Clear versioning**: Semver with documented pre-1.0 expectations
2. **Changelog discipline**: Every PR updates CHANGELOG.md
3. **Release automation**: Tag-triggered publishing with validation
4. **Compatibility guarantees**: Deprecation policy, MSRV commitment
5. **Pre-release support**: Alpha/beta/RC versions for testing

Following this process ensures users can confidently depend on clifford with predictable upgrade paths.
