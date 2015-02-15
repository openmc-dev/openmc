dnl Implement a foreach macro
dnl
divert(-1)
# foreach(x, (item_1, item_2, ..., item_n), stmt)
define(`m4_foreach', `pushdef(`$1', `')_foreach(`$1', `$2', `$3')popdef(`$1')')
define(`_arg1', `$1')
define(`_foreach',
        `ifelse(`$2', `()', ,
                `define(`$1', _arg1$2)$3`'_foreach(`$1', (shift$2), `$3')')')
# traceon(`define', `foreach', `_foreach', `ifelse')
divert`'dnl
