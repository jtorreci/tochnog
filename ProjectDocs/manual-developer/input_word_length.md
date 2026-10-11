# input_word_length

## Where

`input.cc`, static helper `input_word( std::istream &in, char str[],
long int in_comment )`, used by every word read of the input parser:
the initialization part (`input()`), `input_read_string()` (main data
words, include files, `start_define`, `start_arithmetic`, `start_if`
blocks), the comment loop of `input_skip_comment()` (`in_comment` = 1),
`input_runtime()` (`.run` file) and `input_gmsh_read()` (skipped section
words, `in_comment` = 1).

## Implementation

- Reads with `in >> std::setw(MCHAR) >> str`: at most `MCHAR`-1
  characters, the rest of a longer word stays in the stream and becomes
  the next word.
- After the read the next character is peeked with
  `in.rdbuf()->sgetc()` (no stream state change, so the `cin.eof()`
  checks of `input_skip_comment()` behave as before). A non-blank next
  character means the word was split.
- A split word is accepted when `in_comment` is set or when the word
  starts a comment (`str[0] == '('`): the comment parenthesis counting
  sees every piece, also a `)` beyond character 99. Otherwise it is an
  input error (`pri` + `exit(TN_EXIT_STATUS)`).
- Returns 0 at the end of the input, like the former `!(cin >> str)`.

Tests: `validation-suite/test-2026/input_long_comment_word1.dat`
(134- and 108-character words in comments, the second one closes the
comment; rc 0, target met) and `input_long_word_error1.dat` (150-character
record name; expected error `Word too long`). Old GNU: rc 0 by layout
luck / segmentation fault (rc 139); ASan stack-buffer-overflow in
`input_skip_comment()` on the Professional `ground19_water_under_dam`.

## Hardcoded parameters / pending refactorings

- `MCHAR` (100, `tochnog.h`) is the word limit; `MCHAR_WORDS` (60) still
  limits define and arithmetic names separately.
- `range.cc` (`-ra ... -ra` lists, `-from ... -to`) still reads two words
  with a plain `cin >> str` into `char str[MCHAR]`: same overflow, not
  converted in F4 T6a (outside the work unit); move `input_word()` to a
  shared header and use it there.
- Reads into `char*` with `setw` rely on the C++17 `operator>>(istream&,
  char*)`; C++20 removed it (the array-reference overload bounds by
  itself). The build uses the g++ default (gnu++17).
