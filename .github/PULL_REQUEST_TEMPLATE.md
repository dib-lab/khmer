
----------------------------------------

> <u>Message from the khmer maintainers</u>:
>
> *Please provide a description of the pull request above, including a
> reasonable level of detail and references to any relevant threads. After
> creating the  pull request, complete and mark the checklist below. If any item
> is not applicable, feel free to remove it or mark it as complete. If you are
> unsure about any item, feel free to ask about it before or during code
> review.*
>
> *When you are ready for the pull request to be reviewed, please post a comment
> with a message such as "Ready for review!"*

- [ ] Is any new functionality in tested? (This can be checked with
      `make clean diff-cover` or the CodeCov report that is automatically
      generated following a successful CI build.)
- [ ] Was a spellchecker run on the source code and documentation after
      changes were made?
- [ ] Have any changes to the command-line interface been explicitly described?
      Only backwards-compatible additions are allowed without a major version
      increment. Changing file formats also requires a major version number
      increment.
- [ ] Have any substantial changes been documented in `CHANGELOG.md`? See
      [keepachangelog](http://keepachangelog.com/) for more details.
- [ ] Do the changes respect streaming I/O? (Are they tested for streaming I/O?)
