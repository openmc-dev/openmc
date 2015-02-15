dnl
dnl given a list (a, b, c) strip off the brackets:
define(`TOHWM4_dummyarglist',`dnl
substr($1,1,decr(decr(len($1))))`'dnl
')dnl
dnl
dnl given a variable name a, declare it as follows:
define(`TOHWM4_dummyargdecl',`dnl
    character(len=*), intent(in), optional :: $1
')dnl
dnl
dnl use an optional character variable:
define(`TOHWM4_dummyarguse',`dnl
    if (present($1)) call xml_addAttribute(xf, "$1", $1)
')dnl
dnl
define(`TOHWM4_interfacename',module procedure `$1'`$2'`$3')dnl
dnl
define(`TOHWM4_interfacelist', `dnl
     TOHWM4_interfacename(`$1',`$2',`Sca')
     TOHWM4_interfacename(`$1',`$2',`ArrSi')
     TOHWM4_interfacename(`$1',`$2',`ArrSh')
     TOHWM4_interfacename(`$1',`$2',`MatSi')
     TOHWM4_interfacename(`$1',`$2',`MatSh')
')dnl
define(`TOHWM4_interfaceshortlist', `dnl
     TOHWM4_interfacename(`$1',`$2',`Sca')
     TOHWM4_interfacename(`$1',`$2',`Arr')
     TOHWM4_interfacename(`$1',`$2',`Mat')
')dnl
