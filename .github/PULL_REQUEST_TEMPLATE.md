- [ ] Is it mergeable?
- [ ] `make test` Did it pass the tests?
- [ ] `make clean diff-cover` If it introduces new functionality in
  `scripts/` is it tested?
- [ ] `make format diff_pylint_report cppcheck doc pydocstyle` Is it well
  formatted?
- [ ] Did it change the command-line interface? Only backwards-compatible
  additions are allowed without a major version increment. Changing file
  formats also requires a major version number increment.
- [ ] For substantial changes or changes to the command-line interface, is it
  documented in `CHANGELOG.md`? See [keepachangelog](http://keepachangelog.com/)
  for more details.
- [ ] Was a spellchecker run on the source code and documentation after
  changes were made?
- [ ] Do the changes respect streaming IO? (Are they
  tested for streaming IO?)
