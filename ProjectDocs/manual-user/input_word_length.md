# input_word_length

## Description

The input file is read word by word (words are separated by blanks,
tabs or line ends). A word may have at most 99 characters (the parser
buffer `MCHAR` = 100 minus the terminating zero).

- Inside a comment `( ... )` a longer word is allowed: it is read in
  pieces and skipped like the rest of the comment. A typical case is a
  long URL in the header comment of an input file.
- Anywhere else (record names, values, define and arithmetic words, the
  `.run` runtime file, include files) a longer word is an input error:
  the run stops with `Word too long, it starts with <first 99
  characters>` and `maximum allowed length is 99`.

Before F4 T6a (BACKLOG C56) a longer word was written past the parser
buffer (memory corruption): the run could work by chance or crash.

## Usage

Nothing to specify; the rule applies to every input file.

## Example

```
( see https://example.org/a/path/that/is/longer/than/ninety/nine/characters/is/still/fine/inside/a/comment.html )
echo -no
number_of_space_dimensions 1
...
```

## Differences with the Professional version

The Professional 02-08-2026 also accepts long words in comments
(`ground19_water_under_dam`); for a long data word it stops with
`Error in data part.` (exit status 0). The GNU message names the word.
