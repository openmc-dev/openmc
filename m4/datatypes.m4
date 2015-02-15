dnl
define(`TOHWM4_declarationtype', `dnl
ifelse(`$1', `CmplxDp', `complex(dp)',
       `$1', `CmplxSp', `complex(sp)',
       `$1', `RealDp', `real(dp)', 
       `$1', `RealSp', `real(sp)', 
       `$1', `Int', `integer', 
       `$1', `Lg', `logical', 
       `$1', `Ch', `character(len=*)')`'dnl
')dnl
define(`TOHWM4_datatype', `dnl
ifelse(`$1', `CmplxDp', `fpx:complex',
       `$1', `CmplxSp', `fpx:complex',
       `$1', `RealDp', `fpx:real', 
       `$1', `RealSp', `fpx:real', 
       `$1', `Int', `xsd:integer', 
       `$1', `Lg', `xsd:boolean', 
       `$1', `Ch', `xsd:string')`'dnl
')dnl
dnl
define(`TOHWM4_types', `(CmplxDp, CmplxSp, RealDp, RealSp, Int, Lg, Ch)')dnl
