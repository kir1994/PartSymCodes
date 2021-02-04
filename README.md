# Partially symmetric codes
Bounds and code generator for partially symmetric codes (check https://arxiv.org/abs/2001.03790)


### Python script to produce curves similar to Fig. 2
Usage: 
```
python bound.py m t
```
Returns pairs (code dimension, target derivative dimension)

### Partially symmetric monomial code generator 
Usage:
```
MonCodesGen m k t r
```
Prints to standard error stream the description of t-symmetric monomial code of length 2^m and dimension k, which is a subcode of Reed-Muller code RM(r,m). 
First line contains `2^m k`, followed by 2^m-k lines of type `1 sym_id`, where `sym_id` is an index of frozen symbol
