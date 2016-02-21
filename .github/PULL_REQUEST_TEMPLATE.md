- [ ] Is it mergeable?
- [ ] Did it pass the tests?
- [ ] If it introduces new functionality in scripts/ is it tested?
  Check for code coverage with `make clean diff-cover`
- [ ] Is it well formatted? Look at `make pep8`, `make diff_pylint_report`,
  `make cppcheck`, and `make doc` output. Use `make format` and manual
  fixing as needed.
- [ ] Did it change the command-line interface? Only additions are allowed
  without a major version increment. Changing file formats also requires a
  major version number increment.
- [ ] Is it documented in the ChangeLog?
  http://en.wikipedia.org/wiki/Changelog#Format
- [ ] Was a spellchecker run on the source code and documentation after
  changes were made?
- [ ] Do the changes respect streaming IO? (Are they
  tested for streaming IO?)
- [ ] Is the Copyright year up to date?
