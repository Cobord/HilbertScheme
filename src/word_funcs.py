from typing import List, Optional
from hilbert_scheme import TwoDPartition

Nat = int

def is_yamanouchi(y_word : List[Nat]) -> bool:
    n = len(y_word)
    letter_frequencies = {x:0 for x in y_word}
    letter_frequencies.update({max(x-1,0):0 for x in y_word})
    for letter in reversed(y_word):
        letter_frequencies[letter] += 1
        if letter>0 and letter_frequencies[letter]>letter_frequencies[letter-1]:
            return False
    return True

def word_content(a_word : List[Nat]) -> TwoDPartition:
    if not is_yamanouchi(a_word):
        raise ValueError("The word must be Yamanouchi in order to have valid content as a partition of it's length")
    letter_frequencies = {x:0 for x in a_word}
    for letter in a_word:
        letter_frequencies[letter] += 1
    return TwoDPartition(list(letter_frequencies.values()))

def crystal_f(a_word : List[Nat], i : Nat) -> Optional[List[Nat]]:
    num_open_parens = 0
    last_non_cancelled_close_paren = None
    for index, letter in enumerate(a_word):
        if letter == i and num_open_parens==0:
            last_non_cancelled_close_paren = index
        elif letter == i:
            num_open_parens -= 1
        elif letter == i+1:
            num_open_parens += 1
        else:
            pass
    if last_non_cancelled_close_paren is None:
        return None
    else:
        new_word = a_word.copy()
        new_word[last_non_cancelled_close_paren] += 1
        return new_word

def crystal_e(a_word : List[Nat], i : Nat) -> Optional[List[Nat]]:
    for index,letter in enumerate(a_word):
        if letter == i+1:
            temp_word = a_word.copy()
            temp_word[index]-=1
            if crystal_f(temp_word,i)==a_word:
                return temp_word
    return None