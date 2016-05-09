- [ ] Is it mergeable?
- [ ] `make test` Did it pass the tests?
- [ ] `make clean diff-cover` If it introduces new functionality in
  `scripts/` is it tested?
- [ ] `make format diff_pylint_report cppcheck doc pydocstyle` Is it well
  formatted?
- [ ] Did it change the command-line interface? Only additions are allowed
  without a major version increment. Changing file formats also requires a
  major version number increment.
- [ ] Is it documented in the `ChangeLog`?
  http://en.wikipedia.org/wiki/Changelog#Format
- [ ] Was a spellchecker run on the source code and documentation after
  changes were made?
- [ ] Do the changes respect streaming IO? (Are they
  tested for streaming IO?)
- [ ] Is the Copyright year up to date?
