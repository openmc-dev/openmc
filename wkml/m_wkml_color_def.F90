module m_wkml_color_def
#ifndef DUMMYLIB
!Standard Fortran 95 allows 19 continuation line for fixed-form and 39 for free-form
  implicit none
  private

  type color
    integer       :: r, g, b
    character(20) :: name
  end type color

  public :: color
  public :: colorarray

  type(color), parameter :: colorArray0(39) = (/ &
    color(255, 250, 250, 'snow                '), &
    color(248, 248, 255, 'ghost white         '), &
    color(248, 248, 255, 'GhostWhite          '), &
    color(245, 245, 245, 'white smoke         '), &
    color(245, 245, 245, 'WhiteSmoke          '), &
    color(220, 220, 220, 'gainsboro           '), &
    color(255, 250, 240, 'floral white        '), &
    color(255, 250, 240, 'FloralWhite         '), &
    color(253, 245, 230, 'old lace            '), &
    color(253, 245, 230, 'OldLace             '), &
    color(250, 240, 230, 'linen               '), &
    color(250, 235, 215, 'antique white       '), &
    color(250, 235, 215, 'AntiqueWhite        '), &
    color(255, 239, 213, 'papaya whip         '), &
    color(255, 239, 213, 'PapayaWhip          '), &
    color(255, 235, 205, 'blanched almond     '), &
    color(255, 235, 205, 'BlanchedAlmond      '), &
    color(255, 228, 196, 'bisque              '), &
    color(255, 218, 185, 'peach puff          '), &
    color(255, 218, 185, 'PeachPuff           '), &
    color(255, 222, 173, 'navajo white        '), &
    color(255, 222, 173, 'NavajoWhite         '), &
    color(255, 228, 181, 'moccasin            '), &
    color(255, 248, 220, 'cornsilk            '), &
    color(255, 255, 240, 'ivory               '), &
    color(255, 250, 205, 'lemon chiffon       '), &
    color(255, 250, 205, 'LemonChiffon        '), &
    color(255, 245, 238, 'seashell            '), &
    color(240, 255, 240, 'honeydew            '), &
    color(245, 255, 250, 'mint cream          '), &
    color(245, 255, 250, 'MintCream           '), &
    color(240, 255, 255, 'azure               '), &
    color(240, 248, 255, 'alice blue          '), &
    color(240, 248, 255, 'AliceBlue           '), &
    color(230, 230, 250, 'lavender            '), &
    color(255, 240, 245, 'lavender blush      '), &
    color(255, 240, 245, 'LavenderBlush       '), &
    color(255, 228, 225, 'misty rose          '), &
    color(255, 228, 225, 'MistyRose           ')/)

  type(color), parameter :: colorArray1(39) = (/ &
    color(255, 255, 255, 'white               '), &
    color(0, 0, 0, 'black               '), &
    color(47, 79, 79, 'dark slate gray     '), &
    color(47, 79, 79, 'DarkSlateGray       '), &
    color(47, 79, 79, 'dark slate grey     '), &
    color(47, 79, 79, 'DarkSlateGrey       '), &
    color(105, 105, 105, 'dim gray            '), &
    color(105, 105, 105, 'DimGray             '), &
    color(105, 105, 105, 'dim grey            '), &
    color(105, 105, 105, 'DimGrey             '), &
    color(112, 128, 144, 'slate gray          '), &
    color(112, 128, 144, 'SlateGray           '), &
    color(112, 128, 144, 'slate grey          '), &
    color(112, 128, 144, 'SlateGrey           '), &
    color(119, 136, 153, 'light slate gray    '), &
    color(119, 136, 153, 'LightSlateGray      '), &
    color(119, 136, 153, 'light slate grey    '), &
    color(119, 136, 153, 'LightSlateGrey      '), &
    color(190, 190, 190, 'gray                '), &
    color(190, 190, 190, 'grey                '), &
    color(211, 211, 211, 'light grey          '), &
    color(211, 211, 211, 'LightGrey           '), &
    color(211, 211, 211, 'light gray          '), &
    color(211, 211, 211, 'LightGray           '), &
    color(25, 25, 112, 'midnight blue       '), &
    color(25, 25, 112, 'MidnightBlue        '), &
    color(0, 0, 128, 'navy                '), &
    color(0, 0, 128, 'navy blue           '), &
    color(0, 0, 128, 'NavyBlue            '), &
    color(100, 149, 237, 'cornflower blue     '), &
    color(100, 149, 237, 'CornflowerBlue      '), &
    color(72, 61, 139, 'dark slate blue     '), &
    color(72, 61, 139, 'DarkSlateBlue       '), &
    color(106, 90, 205, 'slate blue          '), &
    color(106, 90, 205, 'SlateBlue           '), &
    color(123, 104, 238, 'medium slate blue   '), &
    color(123, 104, 238, 'MediumSlateBlue     '), &
    color(132, 112, 255, 'light slate blue    '), &
    color(132, 112, 255, 'LightSlateBlue      ')/)

  type(color), parameter :: colorArray2(39) = (/ &
    color(0, 0, 205, 'medium blue         '), &
    color(0, 0, 205, 'MediumBlue          '), &
    color(65, 105, 225, 'royal blue          '), &
    color(65, 105, 225, 'RoyalBlue           '), &
    color(0, 0, 255, 'blue                '), &
    color(30, 144, 255, 'dodger blue         '), &
    color(30, 144, 255, 'DodgerBlue          '), &
    color(0, 191, 255, 'deep sky blue       '), &
    color(0, 191, 255, 'DeepSkyBlue         '), &
    color(135, 206, 235, 'sky blue            '), &
    color(135, 206, 235, 'SkyBlue             '), &
    color(135, 206, 250, 'light sky blue      '), &
    color(135, 206, 250, 'LightSkyBlue        '), &
    color(70, 130, 180, 'steel blue          '), &
    color(70, 130, 180, 'SteelBlue           '), &
    color(176, 196, 222, 'light steel blue    '), &
    color(176, 196, 222, 'LightSteelBlue      '), &
    color(173, 216, 230, 'light blue          '), &
    color(173, 216, 230, 'LightBlue           '), &
    color(176, 224, 230, 'powder blue         '), &
    color(176, 224, 230, 'PowderBlue          '), &
    color(175, 238, 238, 'pale turquoise      '), &
    color(175, 238, 238, 'PaleTurquoise       '), &
    color(0, 206, 209, 'dark turquoise      '), &
    color(0, 206, 209, 'DarkTurquoise       '), &
    color(72, 209, 204, 'medium turquoise    '), &
    color(72, 209, 204, 'MediumTurquoise     '), &
    color(64, 224, 208, 'turquoise           '), &
    color(0, 255, 255, 'cyan                '), &
    color(224, 255, 255, 'light cyan          '), &
    color(224, 255, 255, 'LightCyan           '), &
    color(95, 158, 160, 'cadet blue          '), &
    color(95, 158, 160, 'CadetBlue           '), &
    color(102, 205, 170, 'medium aquamarine   '), &
    color(102, 205, 170, 'MediumAquamarine    '), &
    color(127, 255, 212, 'aquamarine          '), &
    color(0, 100, 0, 'dark green          '), &
    color(0, 100, 0, 'DarkGreen           '), &
    color(85, 107, 47, 'dark olive green    ')/)

  type(color), parameter :: colorArray3(39) = (/ &
    color(85, 107, 47, 'DarkOliveGreen      '), &
    color(143, 188, 143, 'dark sea green      '), &
    color(143, 188, 143, 'DarkSeaGreen        '), &
    color(46, 139, 87, 'sea green           '), &
    color(46, 139, 87, 'SeaGreen            '), &
    color(60, 179, 113, 'medium sea green    '), &
    color(60, 179, 113, 'MediumSeaGreen      '), &
    color(32, 178, 170, 'light sea green     '), &
    color(32, 178, 170, 'LightSeaGreen       '), &
    color(152, 251, 152, 'pale green          '), &
    color(152, 251, 152, 'PaleGreen           '), &
    color(0, 255, 127, 'spring green        '), &
    color(0, 255, 127, 'SpringGreen         '), &
    color(124, 252, 0, 'lawn green          '), &
    color(124, 252, 0, 'LawnGreen           '), &
    color(0, 255, 0, 'green               '), &
    color(127, 255, 0, 'chartreuse          '), &
    color(0, 250, 154, 'medium spring green '), &
    color(0, 250, 154, 'MediumSpringGreen   '), &
    color(173, 255, 47, 'green yellow        '), &
    color(173, 255, 47, 'GreenYellow         '), &
    color(50, 205, 50, 'lime green          '), &
    color(50, 205, 50, 'LimeGreen           '), &
    color(154, 205, 50, 'yellow green        '), &
    color(154, 205, 50, 'YellowGreen         '), &
    color(34, 139, 34, 'forest green        '), &
    color(34, 139, 34, 'ForestGreen         '), &
    color(107, 142, 35, 'olive drab          '), &
    color(107, 142, 35, 'OliveDrab           '), &
    color(189, 183, 107, 'dark khaki          '), &
    color(189, 183, 107, 'DarkKhaki           '), &
    color(240, 230, 140, 'khaki               '), &
    color(238, 232, 170, 'pale goldenrod      '), &
    color(238, 232, 170, 'PaleGoldenrod       '), &
    color(250, 250, 210, 'LightGoldenrodYellow'), &
    color(255, 255, 224, 'light yellow        '), &
    color(255, 255, 224, 'LightYellow         '), &
    color(255, 255, 0, 'yellow              '), &
    color(255, 215, 0, 'gold                ')/)

  type(color), parameter :: colorArray4(39) = (/ &
    color(238, 221, 130, 'light goldenrod     '), &
    color(238, 221, 130, 'LightGoldenrod      '), &
    color(218, 165, 32, 'goldenrod           '), &
    color(184, 134, 11, 'dark goldenrod      '), &
    color(184, 134, 11, 'DarkGoldenrod       '), &
    color(188, 143, 143, 'rosy brown          '), &
    color(188, 143, 143, 'RosyBrown           '), &
    color(205, 92, 92, 'indian red          '), &
    color(205, 92, 92, 'IndianRed           '), &
    color(139, 69, 19, 'saddle brown        '), &
    color(139, 69, 19, 'SaddleBrown         '), &
    color(160, 82, 45, 'sienna              '), &
    color(205, 133, 63, 'peru                '), &
    color(222, 184, 135, 'burlywood           '), &
    color(245, 245, 220, 'beige               '), &
    color(245, 222, 179, 'wheat               '), &
    color(244, 164, 96, 'sandy brown         '), &
    color(244, 164, 96, 'SandyBrown          '), &
    color(210, 180, 140, 'tan                 '), &
    color(210, 105, 30, 'chocolate           '), &
    color(178, 34, 34, 'firebrick           '), &
    color(165, 42, 42, 'brown               '), &
    color(233, 150, 122, 'dark salmon         '), &
    color(233, 150, 122, 'DarkSalmon          '), &
    color(250, 128, 114, 'salmon              '), &
    color(255, 160, 122, 'light salmon        '), &
    color(255, 160, 122, 'LightSalmon         '), &
    color(255, 165, 0, 'orange              '), &
    color(255, 140, 0, 'dark orange         '), &
    color(255, 140, 0, 'DarkOrange          '), &
    color(255, 127, 80, 'coral               '), &
    color(240, 128, 128, 'light coral         '), &
    color(240, 128, 128, 'LightCoral          '), &
    color(255, 99, 71, 'tomato              '), &
    color(255, 69, 0, 'orange red          '), &
    color(255, 69, 0, 'OrangeRed           '), &
    color(255, 0, 0, 'red                 '), &
    color(255, 105, 180, 'hot pink            '), &
    color(255, 105, 180, 'HotPink             ')/)

  type(color), parameter :: colorArray5(39) = (/ &
    color(255, 20, 147, 'deep pink           '), &
    color(255, 20, 147, 'DeepPink            '), &
    color(255, 192, 203, 'pink                '), &
    color(255, 182, 193, 'light pink          '), &
    color(255, 182, 193, 'LightPink           '), &
    color(219, 112, 147, 'pale violet red     '), &
    color(219, 112, 147, 'PaleVioletRed       '), &
    color(176, 48, 96, 'maroon              '), &
    color(199, 21, 133, 'medium violet red   '), &
    color(199, 21, 133, 'MediumVioletRed     '), &
    color(208, 32, 144, 'violet red          '), &
    color(208, 32, 144, 'VioletRed           '), &
    color(255, 0, 255, 'magenta             '), &
    color(238, 130, 238, 'violet              '), &
    color(221, 160, 221, 'plum                '), &
    color(218, 112, 214, 'orchid              '), &
    color(186, 85, 211, 'medium orchid       '), &
    color(186, 85, 211, 'MediumOrchid        '), &
    color(153, 50, 204, 'dark orchid         '), &
    color(153, 50, 204, 'DarkOrchid          '), &
    color(148, 0, 211, 'dark violet         '), &
    color(148, 0, 211, 'DarkViolet          '), &
    color(138, 43, 226, 'blue violet         '), &
    color(138, 43, 226, 'BlueViolet          '), &
    color(160, 32, 240, 'purple              '), &
    color(147, 112, 219, 'medium purple       '), &
    color(147, 112, 219, 'MediumPurple        '), &
    color(216, 191, 216, 'thistle             '), &
    color(255, 250, 250, 'snow1               '), &
    color(238, 233, 233, 'snow2               '), &
    color(205, 201, 201, 'snow3               '), &
    color(139, 137, 137, 'snow4               '), &
    color(255, 245, 238, 'seashell1           '), &
    color(238, 229, 222, 'seashell2           '), &
    color(205, 197, 191, 'seashell3           '), &
    color(139, 134, 130, 'seashell4           '), &
    color(255, 239, 219, 'AntiqueWhite1       '), &
    color(238, 223, 204, 'AntiqueWhite2       '), &
    color(205, 192, 176, 'AntiqueWhite3       ')/)

  type(color), parameter :: colorArray6(39) = (/ &
    color(139, 131, 120, 'AntiqueWhite4       '), &
    color(255, 228, 196, 'bisque1             '), &
    color(238, 213, 183, 'bisque2             '), &
    color(205, 183, 158, 'bisque3             '), &
    color(139, 125, 107, 'bisque4             '), &
    color(255, 218, 185, 'PeachPuff1          '), &
    color(238, 203, 173, 'PeachPuff2          '), &
    color(205, 175, 149, 'PeachPuff3          '), &
    color(139, 119, 101, 'PeachPuff4          '), &
    color(255, 222, 173, 'NavajoWhite1        '), &
    color(238, 207, 161, 'NavajoWhite2        '), &
    color(205, 179, 139, 'NavajoWhite3        '), &
    color(139, 121, 94, 'NavajoWhite4        '), &
    color(255, 250, 205, 'LemonChiffon1       '), &
    color(238, 233, 191, 'LemonChiffon2       '), &
    color(205, 201, 165, 'LemonChiffon3       '), &
    color(139, 137, 112, 'LemonChiffon4       '), &
    color(255, 248, 220, 'cornsilk1           '), &
    color(238, 232, 205, 'cornsilk2           '), &
    color(205, 200, 177, 'cornsilk3           '), &
    color(139, 136, 120, 'cornsilk4           '), &
    color(255, 255, 240, 'ivory1              '), &
    color(238, 238, 224, 'ivory2              '), &
    color(205, 205, 193, 'ivory3              '), &
    color(139, 139, 131, 'ivory4              '), &
    color(240, 255, 240, 'honeydew1           '), &
    color(224, 238, 224, 'honeydew2           '), &
    color(193, 205, 193, 'honeydew3           '), &
    color(131, 139, 131, 'honeydew4           '), &
    color(255, 240, 245, 'LavenderBlush1      '), &
    color(238, 224, 229, 'LavenderBlush2      '), &
    color(205, 193, 197, 'LavenderBlush3      '), &
    color(139, 131, 134, 'LavenderBlush4      '), &
    color(255, 228, 225, 'MistyRose1          '), &
    color(238, 213, 210, 'MistyRose2          '), &
    color(205, 183, 181, 'MistyRose3          '), &
    color(139, 125, 123, 'MistyRose4          '), &
    color(240, 255, 255, 'azure1              '), &
    color(224, 238, 238, 'azure2              ')/)

  type(color), parameter :: colorArray7(39) = (/ &
    color(193, 205, 205, 'azure3              '), &
    color(131, 139, 139, 'azure4              '), &
    color(131, 111, 255, 'SlateBlue1          '), &
    color(122, 103, 238, 'SlateBlue2          '), &
    color(105, 89, 205, 'SlateBlue3          '), &
    color(71, 60, 139, 'SlateBlue4          '), &
    color(72, 118, 255, 'RoyalBlue1          '), &
    color(67, 110, 238, 'RoyalBlue2          '), &
    color(58, 95, 205, 'RoyalBlue3          '), &
    color(39, 64, 139, 'RoyalBlue4          '), &
    color(0, 0, 255, 'blue1               '), &
    color(0, 0, 238, 'blue2               '), &
    color(0, 0, 205, 'blue3               '), &
    color(0, 0, 139, 'blue4               '), &
    color(30, 144, 255, 'DodgerBlue1         '), &
    color(28, 134, 238, 'DodgerBlue2         '), &
    color(24, 116, 205, 'DodgerBlue3         '), &
    color(16, 78, 139, 'DodgerBlue4         '), &
    color(99, 184, 255, 'SteelBlue1          '), &
    color(92, 172, 238, 'SteelBlue2          '), &
    color(79, 148, 205, 'SteelBlue3          '), &
    color(54, 100, 139, 'SteelBlue4          '), &
    color(0, 191, 255, 'DeepSkyBlue1        '), &
    color(0, 178, 238, 'DeepSkyBlue2        '), &
    color(0, 154, 205, 'DeepSkyBlue3        '), &
    color(0, 104, 139, 'DeepSkyBlue4        '), &
    color(135, 206, 255, 'SkyBlue1            '), &
    color(126, 192, 238, 'SkyBlue2            '), &
    color(108, 166, 205, 'SkyBlue3            '), &
    color(74, 112, 139, 'SkyBlue4            '), &
    color(176, 226, 255, 'LightSkyBlue1       '), &
    color(164, 211, 238, 'LightSkyBlue2       '), &
    color(141, 182, 205, 'LightSkyBlue3       '), &
    color(96, 123, 139, 'LightSkyBlue4       '), &
    color(198, 226, 255, 'SlateGray1          '), &
    color(185, 211, 238, 'SlateGray2          '), &
    color(159, 182, 205, 'SlateGray3          '), &
    color(108, 123, 139, 'SlateGray4          '), &
    color(202, 225, 255, 'LightSteelBlue1     ')/)

  type(color), parameter :: colorArray8(39) = (/ &
    color(188, 210, 238, 'LightSteelBlue2     '), &
    color(162, 181, 205, 'LightSteelBlue3     '), &
    color(110, 123, 139, 'LightSteelBlue4     '), &
    color(191, 239, 255, 'LightBlue1          '), &
    color(178, 223, 238, 'LightBlue2          '), &
    color(154, 192, 205, 'LightBlue3          '), &
    color(104, 131, 139, 'LightBlue4          '), &
    color(224, 255, 255, 'LightCyan1          '), &
    color(209, 238, 238, 'LightCyan2          '), &
    color(180, 205, 205, 'LightCyan3          '), &
    color(122, 139, 139, 'LightCyan4          '), &
    color(187, 255, 255, 'PaleTurquoise1      '), &
    color(174, 238, 238, 'PaleTurquoise2      '), &
    color(150, 205, 205, 'PaleTurquoise3      '), &
    color(102, 139, 139, 'PaleTurquoise4      '), &
    color(152, 245, 255, 'CadetBlue1          '), &
    color(142, 229, 238, 'CadetBlue2          '), &
    color(122, 197, 205, 'CadetBlue3          '), &
    color(83, 134, 139, 'CadetBlue4          '), &
    color(0, 245, 255, 'turquoise1          '), &
    color(0, 229, 238, 'turquoise2          '), &
    color(0, 197, 205, 'turquoise3          '), &
    color(0, 134, 139, 'turquoise4          '), &
    color(0, 255, 255, 'cyan1               '), &
    color(0, 238, 238, 'cyan2               '), &
    color(0, 205, 205, 'cyan3               '), &
    color(0, 139, 139, 'cyan4               '), &
    color(151, 255, 255, 'DarkSlateGray1      '), &
    color(141, 238, 238, 'DarkSlateGray2      '), &
    color(121, 205, 205, 'DarkSlateGray3      '), &
    color(82, 139, 139, 'DarkSlateGray4      '), &
    color(127, 255, 212, 'aquamarine1         '), &
    color(118, 238, 198, 'aquamarine2         '), &
    color(102, 205, 170, 'aquamarine3         '), &
    color(69, 139, 116, 'aquamarine4         '), &
    color(193, 255, 193, 'DarkSeaGreen1       '), &
    color(180, 238, 180, 'DarkSeaGreen2       '), &
    color(155, 205, 155, 'DarkSeaGreen3       '), &
    color(105, 139, 105, 'DarkSeaGreen4       ')/)

  type(color), parameter :: colorArray9(39) = (/ &
    color(84, 255, 159, 'SeaGreen1           '), &
    color(78, 238, 148, 'SeaGreen2           '), &
    color(67, 205, 128, 'SeaGreen3           '), &
    color(46, 139, 87, 'SeaGreen4           '), &
    color(154, 255, 154, 'PaleGreen1          '), &
    color(144, 238, 144, 'PaleGreen2          '), &
    color(124, 205, 124, 'PaleGreen3          '), &
    color(84, 139, 84, 'PaleGreen4          '), &
    color(0, 255, 127, 'SpringGreen1        '), &
    color(0, 238, 118, 'SpringGreen2        '), &
    color(0, 205, 102, 'SpringGreen3        '), &
    color(0, 139, 69, 'SpringGreen4        '), &
    color(0, 255, 0, 'green1              '), &
    color(0, 238, 0, 'green2              '), &
    color(0, 205, 0, 'green3              '), &
    color(0, 139, 0, 'green4              '), &
    color(127, 255, 0, 'chartreuse1         '), &
    color(118, 238, 0, 'chartreuse2         '), &
    color(102, 205, 0, 'chartreuse3         '), &
    color(69, 139, 0, 'chartreuse4         '), &
    color(192, 255, 62, 'OliveDrab1          '), &
    color(179, 238, 58, 'OliveDrab2          '), &
    color(154, 205, 50, 'OliveDrab3          '), &
    color(105, 139, 34, 'OliveDrab4          '), &
    color(202, 255, 112, 'DarkOliveGreen1     '), &
    color(188, 238, 104, 'DarkOliveGreen2     '), &
    color(162, 205, 90, 'DarkOliveGreen3     '), &
    color(110, 139, 61, 'DarkOliveGreen4     '), &
    color(255, 246, 143, 'khaki1              '), &
    color(238, 230, 133, 'khaki2              '), &
    color(205, 198, 115, 'khaki3              '), &
    color(139, 134, 78, 'khaki4              '), &
    color(255, 236, 139, 'LightGoldenrod1     '), &
    color(238, 220, 130, 'LightGoldenrod2     '), &
    color(205, 190, 112, 'LightGoldenrod3     '), &
    color(139, 129, 76, 'LightGoldenrod4     '), &
    color(255, 255, 224, 'LightYellow1        '), &
    color(238, 238, 209, 'LightYellow2        '), &
    color(205, 205, 180, 'LightYellow3        ')/)

  type(color), parameter :: colorArray10(39) = (/ &
    color(139, 139, 122, 'LightYellow4        '), &
    color(255, 255, 0, 'yellow1             '), &
    color(238, 238, 0, 'yellow2             '), &
    color(205, 205, 0, 'yellow3             '), &
    color(139, 139, 0, 'yellow4             '), &
    color(255, 215, 0, 'gold1               '), &
    color(238, 201, 0, 'gold2               '), &
    color(205, 173, 0, 'gold3               '), &
    color(139, 117, 0, 'gold4               '), &
    color(255, 193, 37, 'goldenrod1          '), &
    color(238, 180, 34, 'goldenrod2          '), &
    color(205, 155, 29, 'goldenrod3          '), &
    color(139, 105, 20, 'goldenrod4          '), &
    color(255, 185, 15, 'DarkGoldenrod1      '), &
    color(238, 173, 14, 'DarkGoldenrod2      '), &
    color(205, 149, 12, 'DarkGoldenrod3      '), &
    color(139, 101, 8, 'DarkGoldenrod4      '), &
    color(255, 193, 193, 'RosyBrown1          '), &
    color(238, 180, 180, 'RosyBrown2          '), &
    color(205, 155, 155, 'RosyBrown3          '), &
    color(139, 105, 105, 'RosyBrown4          '), &
    color(255, 106, 106, 'IndianRed1          '), &
    color(238, 99, 99, 'IndianRed2          '), &
    color(205, 85, 85, 'IndianRed3          '), &
    color(139, 58, 58, 'IndianRed4          '), &
    color(255, 130, 71, 'sienna1             '), &
    color(238, 121, 66, 'sienna2             '), &
    color(205, 104, 57, 'sienna3             '), &
    color(139, 71, 38, 'sienna4             '), &
    color(255, 211, 155, 'burlywood1          '), &
    color(238, 197, 145, 'burlywood2          '), &
    color(205, 170, 125, 'burlywood3          '), &
    color(139, 115, 85, 'burlywood4          '), &
    color(255, 231, 186, 'wheat1              '), &
    color(238, 216, 174, 'wheat2              '), &
    color(205, 186, 150, 'wheat3              '), &
    color(139, 126, 102, 'wheat4              '), &
    color(255, 165, 79, 'tan1                '), &
    color(238, 154, 73, 'tan2                ')/)

  type(color), parameter :: colorArray11(39) = (/ &
    color(205, 133, 63, 'tan3                '), &
    color(139, 90, 43, 'tan4                '), &
    color(255, 127, 36, 'chocolate1          '), &
    color(238, 118, 33, 'chocolate2          '), &
    color(205, 102, 29, 'chocolate3          '), &
    color(139, 69, 19, 'chocolate4          '), &
    color(255, 48, 48, 'firebrick1          '), &
    color(238, 44, 44, 'firebrick2          '), &
    color(205, 38, 38, 'firebrick3          '), &
    color(139, 26, 26, 'firebrick4          '), &
    color(255, 64, 64, 'brown1              '), &
    color(238, 59, 59, 'brown2              '), &
    color(205, 51, 51, 'brown3              '), &
    color(139, 35, 35, 'brown4              '), &
    color(255, 140, 105, 'salmon1             '), &
    color(238, 130, 98, 'salmon2             '), &
    color(205, 112, 84, 'salmon3             '), &
    color(139, 76, 57, 'salmon4             '), &
    color(255, 160, 122, 'LightSalmon1        '), &
    color(238, 149, 114, 'LightSalmon2        '), &
    color(205, 129, 98, 'LightSalmon3        '), &
    color(139, 87, 66, 'LightSalmon4        '), &
    color(255, 165, 0, 'orange1             '), &
    color(238, 154, 0, 'orange2             '), &
    color(205, 133, 0, 'orange3             '), &
    color(139, 90, 0, 'orange4             '), &
    color(255, 127, 0, 'DarkOrange1         '), &
    color(238, 118, 0, 'DarkOrange2         '), &
    color(205, 102, 0, 'DarkOrange3         '), &
    color(139, 69, 0, 'DarkOrange4         '), &
    color(255, 114, 86, 'coral1              '), &
    color(238, 106, 80, 'coral2              '), &
    color(205, 91, 69, 'coral3              '), &
    color(139, 62, 47, 'coral4              '), &
    color(255, 99, 71, 'tomato1             '), &
    color(238, 92, 66, 'tomato2             '), &
    color(205, 79, 57, 'tomato3             '), &
    color(139, 54, 38, 'tomato4             '), &
    color(255, 69, 0, 'OrangeRed1          ')/)

  type(color), parameter :: colorArray12(39) = (/ &
    color(238, 64, 0, 'OrangeRed2          '), &
    color(205, 55, 0, 'OrangeRed3          '), &
    color(139, 37, 0, 'OrangeRed4          '), &
    color(255, 0, 0, 'red1                '), &
    color(238, 0, 0, 'red2                '), &
    color(205, 0, 0, 'red3                '), &
    color(139, 0, 0, 'red4                '), &
    color(255, 20, 147, 'DeepPink1           '), &
    color(238, 18, 137, 'DeepPink2           '), &
    color(205, 16, 118, 'DeepPink3           '), &
    color(139, 10, 80, 'DeepPink4           '), &
    color(255, 110, 180, 'HotPink1            '), &
    color(238, 106, 167, 'HotPink2            '), &
    color(205, 96, 144, 'HotPink3            '), &
    color(139, 58, 98, 'HotPink4            '), &
    color(255, 181, 197, 'pink1               '), &
    color(238, 169, 184, 'pink2               '), &
    color(205, 145, 158, 'pink3               '), &
    color(139, 99, 108, 'pink4               '), &
    color(255, 174, 185, 'LightPink1          '), &
    color(238, 162, 173, 'LightPink2          '), &
    color(205, 140, 149, 'LightPink3          '), &
    color(139, 95, 101, 'LightPink4          '), &
    color(255, 130, 171, 'PaleVioletRed1      '), &
    color(238, 121, 159, 'PaleVioletRed2      '), &
    color(205, 104, 137, 'PaleVioletRed3      '), &
    color(139, 71, 93, 'PaleVioletRed4      '), &
    color(255, 52, 179, 'maroon1             '), &
    color(238, 48, 167, 'maroon2             '), &
    color(205, 41, 144, 'maroon3             '), &
    color(139, 28, 98, 'maroon4             '), &
    color(255, 62, 150, 'VioletRed1          '), &
    color(238, 58, 140, 'VioletRed2          '), &
    color(205, 50, 120, 'VioletRed3          '), &
    color(139, 34, 82, 'VioletRed4          '), &
    color(255, 0, 255, 'magenta1            '), &
    color(238, 0, 238, 'magenta2            '), &
    color(205, 0, 205, 'magenta3            '), &
    color(139, 0, 139, 'magenta4            ')/)

  type(color), parameter :: colorArray13(39) = (/ &
    color(255, 131, 250, 'orchid1             '), &
    color(238, 122, 233, 'orchid2             '), &
    color(205, 105, 201, 'orchid3             '), &
    color(139, 71, 137, 'orchid4             '), &
    color(255, 187, 255, 'plum1               '), &
    color(238, 174, 238, 'plum2               '), &
    color(205, 150, 205, 'plum3               '), &
    color(139, 102, 139, 'plum4               '), &
    color(224, 102, 255, 'MediumOrchid1       '), &
    color(209, 95, 238, 'MediumOrchid2       '), &
    color(180, 82, 205, 'MediumOrchid3       '), &
    color(122, 55, 139, 'MediumOrchid4       '), &
    color(191, 62, 255, 'DarkOrchid1         '), &
    color(178, 58, 238, 'DarkOrchid2         '), &
    color(154, 50, 205, 'DarkOrchid3         '), &
    color(104, 34, 139, 'DarkOrchid4         '), &
    color(155, 48, 255, 'purple1             '), &
    color(145, 44, 238, 'purple2             '), &
    color(125, 38, 205, 'purple3             '), &
    color(85, 26, 139, 'purple4             '), &
    color(171, 130, 255, 'MediumPurple1       '), &
    color(159, 121, 238, 'MediumPurple2       '), &
    color(137, 104, 205, 'MediumPurple3       '), &
    color(93, 71, 139, 'MediumPurple4       '), &
    color(255, 225, 255, 'thistle1            '), &
    color(238, 210, 238, 'thistle2            '), &
    color(205, 181, 205, 'thistle3            '), &
    color(139, 123, 139, 'thistle4            '), &
    color(0, 0, 0, 'gray0               '), &
    color(0, 0, 0, 'grey0               '), &
    color(3, 3, 3, 'gray1               '), &
    color(3, 3, 3, 'grey1               '), &
    color(5, 5, 5, 'gray2               '), &
    color(5, 5, 5, 'grey2               '), &
    color(8, 8, 8, 'gray3               '), &
    color(8, 8, 8, 'grey3               '), &
    color(10, 10, 10, 'gray4               '), &
    color(10, 10, 10, 'grey4               '), &
    color(13, 13, 13, 'gray5               ')/)

  type(color), parameter :: colorArray14(39) = (/ &
    color(13, 13, 13, 'grey5               '), &
    color(15, 15, 15, 'gray6               '), &
    color(15, 15, 15, 'grey6               '), &
    color(18, 18, 18, 'gray7               '), &
    color(18, 18, 18, 'grey7               '), &
    color(20, 20, 20, 'gray8               '), &
    color(20, 20, 20, 'grey8               '), &
    color(23, 23, 23, 'gray9               '), &
    color(23, 23, 23, 'grey9               '), &
    color(26, 26, 26, 'gray10              '), &
    color(26, 26, 26, 'grey10              '), &
    color(28, 28, 28, 'gray11              '), &
    color(28, 28, 28, 'grey11              '), &
    color(31, 31, 31, 'gray12              '), &
    color(31, 31, 31, 'grey12              '), &
    color(33, 33, 33, 'gray13              '), &
    color(33, 33, 33, 'grey13              '), &
    color(36, 36, 36, 'gray14              '), &
    color(36, 36, 36, 'grey14              '), &
    color(38, 38, 38, 'gray15              '), &
    color(38, 38, 38, 'grey15              '), &
    color(41, 41, 41, 'gray16              '), &
    color(41, 41, 41, 'grey16              '), &
    color(43, 43, 43, 'gray17              '), &
    color(43, 43, 43, 'grey17              '), &
    color(46, 46, 46, 'gray18              '), &
    color(46, 46, 46, 'grey18              '), &
    color(48, 48, 48, 'gray19              '), &
    color(48, 48, 48, 'grey19              '), &
    color(51, 51, 51, 'gray20              '), &
    color(51, 51, 51, 'grey20              '), &
    color(54, 54, 54, 'gray21              '), &
    color(54, 54, 54, 'grey21              '), &
    color(56, 56, 56, 'gray22              '), &
    color(56, 56, 56, 'grey22              '), &
    color(59, 59, 59, 'gray23              '), &
    color(59, 59, 59, 'grey23              '), &
    color(61, 61, 61, 'gray24              '), &
    color(61, 61, 61, 'grey24              ')/)

  type(color), parameter :: colorArray15(39) = (/ &
    color(64, 64, 64, 'gray25              '), &
    color(64, 64, 64, 'grey25              '), &
    color(66, 66, 66, 'gray26              '), &
    color(66, 66, 66, 'grey26              '), &
    color(69, 69, 69, 'gray27              '), &
    color(69, 69, 69, 'grey27              '), &
    color(71, 71, 71, 'gray28              '), &
    color(71, 71, 71, 'grey28              '), &
    color(74, 74, 74, 'gray29              '), &
    color(74, 74, 74, 'grey29              '), &
    color(77, 77, 77, 'gray30              '), &
    color(77, 77, 77, 'grey30              '), &
    color(79, 79, 79, 'gray31              '), &
    color(79, 79, 79, 'grey31              '), &
    color(82, 82, 82, 'gray32              '), &
    color(82, 82, 82, 'grey32              '), &
    color(84, 84, 84, 'gray33              '), &
    color(84, 84, 84, 'grey33              '), &
    color(87, 87, 87, 'gray34              '), &
    color(87, 87, 87, 'grey34              '), &
    color(89, 89, 89, 'gray35              '), &
    color(89, 89, 89, 'grey35              '), &
    color(92, 92, 92, 'gray36              '), &
    color(92, 92, 92, 'grey36              '), &
    color(94, 94, 94, 'gray37              '), &
    color(94, 94, 94, 'grey37              '), &
    color(97, 97, 97, 'gray38              '), &
    color(97, 97, 97, 'grey38              '), &
    color(99, 99, 99, 'gray39              '), &
    color(99, 99, 99, 'grey39              '), &
    color(102, 102, 102, 'gray40              '), &
    color(102, 102, 102, 'grey40              '), &
    color(105, 105, 105, 'gray41              '), &
    color(105, 105, 105, 'grey41              '), &
    color(107, 107, 107, 'gray42              '), &
    color(107, 107, 107, 'grey42              '), &
    color(110, 110, 110, 'gray43              '), &
    color(110, 110, 110, 'grey43              '), &
    color(112, 112, 112, 'gray44              ')/)

  type(color), parameter :: colorArray16(39) = (/ &
    color(112, 112, 112, 'grey44              '), &
    color(115, 115, 115, 'gray45              '), &
    color(115, 115, 115, 'grey45              '), &
    color(117, 117, 117, 'gray46              '), &
    color(117, 117, 117, 'grey46              '), &
    color(120, 120, 120, 'gray47              '), &
    color(120, 120, 120, 'grey47              '), &
    color(122, 122, 122, 'gray48              '), &
    color(122, 122, 122, 'grey48              '), &
    color(125, 125, 125, 'gray49              '), &
    color(125, 125, 125, 'grey49              '), &
    color(127, 127, 127, 'gray50              '), &
    color(127, 127, 127, 'grey50              '), &
    color(130, 130, 130, 'gray51              '), &
    color(130, 130, 130, 'grey51              '), &
    color(133, 133, 133, 'gray52              '), &
    color(133, 133, 133, 'grey52              '), &
    color(135, 135, 135, 'gray53              '), &
    color(135, 135, 135, 'grey53              '), &
    color(138, 138, 138, 'gray54              '), &
    color(138, 138, 138, 'grey54              '), &
    color(140, 140, 140, 'gray55              '), &
    color(140, 140, 140, 'grey55              '), &
    color(143, 143, 143, 'gray56              '), &
    color(143, 143, 143, 'grey56              '), &
    color(145, 145, 145, 'gray57              '), &
    color(145, 145, 145, 'grey57              '), &
    color(148, 148, 148, 'gray58              '), &
    color(148, 148, 148, 'grey58              '), &
    color(150, 150, 150, 'gray59              '), &
    color(150, 150, 150, 'grey59              '), &
    color(153, 153, 153, 'gray60              '), &
    color(153, 153, 153, 'grey60              '), &
    color(156, 156, 156, 'gray61              '), &
    color(156, 156, 156, 'grey61              '), &
    color(158, 158, 158, 'gray62              '), &
    color(158, 158, 158, 'grey62              '), &
    color(161, 161, 161, 'gray63              '), &
    color(161, 161, 161, 'grey63              ')/)

  type(color), parameter :: colorArray17(39) = (/ &
    color(163, 163, 163, 'gray64              '), &
    color(163, 163, 163, 'grey64              '), &
    color(166, 166, 166, 'gray65              '), &
    color(166, 166, 166, 'grey65              '), &
    color(168, 168, 168, 'gray66              '), &
    color(168, 168, 168, 'grey66              '), &
    color(171, 171, 171, 'gray67              '), &
    color(171, 171, 171, 'grey67              '), &
    color(173, 173, 173, 'gray68              '), &
    color(173, 173, 173, 'grey68              '), &
    color(176, 176, 176, 'gray69              '), &
    color(176, 176, 176, 'grey69              '), &
    color(179, 179, 179, 'gray70              '), &
    color(179, 179, 179, 'grey70              '), &
    color(181, 181, 181, 'gray71              '), &
    color(181, 181, 181, 'grey71              '), &
    color(184, 184, 184, 'gray72              '), &
    color(184, 184, 184, 'grey72              '), &
    color(186, 186, 186, 'gray73              '), &
    color(186, 186, 186, 'grey73              '), &
    color(189, 189, 189, 'gray74              '), &
    color(189, 189, 189, 'grey74              '), &
    color(191, 191, 191, 'gray75              '), &
    color(191, 191, 191, 'grey75              '), &
    color(194, 194, 194, 'gray76              '), &
    color(194, 194, 194, 'grey76              '), &
    color(196, 196, 196, 'gray77              '), &
    color(196, 196, 196, 'grey77              '), &
    color(199, 199, 199, 'gray78              '), &
    color(199, 199, 199, 'grey78              '), &
    color(201, 201, 201, 'gray79              '), &
    color(201, 201, 201, 'grey79              '), &
    color(204, 204, 204, 'gray80              '), &
    color(204, 204, 204, 'grey80              '), &
    color(207, 207, 207, 'gray81              '), &
    color(207, 207, 207, 'grey81              '), &
    color(209, 209, 209, 'gray82              '), &
    color(209, 209, 209, 'grey82              '), &
    color(212, 212, 212, 'gray83              ')/)

  type(color), parameter :: colorArray18(39) = (/ &
    color(212, 212, 212, 'grey83              '), &
    color(214, 214, 214, 'gray84              '), &
    color(214, 214, 214, 'grey84              '), &
    color(217, 217, 217, 'gray85              '), &
    color(217, 217, 217, 'grey85              '), &
    color(219, 219, 219, 'gray86              '), &
    color(219, 219, 219, 'grey86              '), &
    color(222, 222, 222, 'gray87              '), &
    color(222, 222, 222, 'grey87              '), &
    color(224, 224, 224, 'gray88              '), &
    color(224, 224, 224, 'grey88              '), &
    color(227, 227, 227, 'gray89              '), &
    color(227, 227, 227, 'grey89              '), &
    color(229, 229, 229, 'gray90              '), &
    color(229, 229, 229, 'grey90              '), &
    color(232, 232, 232, 'gray91              '), &
    color(232, 232, 232, 'grey91              '), &
    color(235, 235, 235, 'gray92              '), &
    color(235, 235, 235, 'grey92              '), &
    color(237, 237, 237, 'gray93              '), &
    color(237, 237, 237, 'grey93              '), &
    color(240, 240, 240, 'gray94              '), &
    color(240, 240, 240, 'grey94              '), &
    color(242, 242, 242, 'gray95              '), &
    color(242, 242, 242, 'grey95              '), &
    color(245, 245, 245, 'gray96              '), &
    color(245, 245, 245, 'grey96              '), &
    color(247, 247, 247, 'gray97              '), &
    color(247, 247, 247, 'grey97              '), &
    color(250, 250, 250, 'gray98              '), &
    color(250, 250, 250, 'grey98              '), &
    color(252, 252, 252, 'gray99              '), &
    color(252, 252, 252, 'grey99              '), &
    color(255, 255, 255, 'gray100             '), &
    color(255, 255, 255, 'grey100             '), &
    color(169, 169, 169, 'dark grey           '), &
    color(169, 169, 169, 'DarkGrey            '), &
    color(169, 169, 169, 'dark gray           '), &
    color(169, 169, 169, 'DarkGray            ')/)


  type(color), parameter :: colorArray19(10) = (/ &
    color(0, 0, 139, 'dark blue           '), &
    color(0, 0, 139, 'DarkBlue            '), &
    color(0, 139, 139, 'dark cyan           '), &
    color(0, 139, 139, 'DarkCyan            '), &
    color(139, 0, 139, 'dark magenta        '), &
    color(139, 0, 139, 'DarkMagenta         '), &
    color(139, 0, 0, 'dark red            '), &
    color(139, 0, 0, 'DarkRed             '), &
    color(144, 238, 144, 'light green         '), &
    color(144, 238, 144, 'LightGreen          ')/)

  type(color), parameter :: colorArray(751) = (/ &
    colorArray0, &
    colorArray1, &
    colorArray2, &
    colorArray3, &
    colorArray4, &
    colorArray5, &
    colorArray6, &
    colorArray7, &
    colorArray8, &
    colorArray9, &
    colorArray10, &
    colorArray11, &
    colorArray12, &
    colorArray13, &
    colorArray14, &
    colorArray15, &
    colorArray16, &
    colorArray17, &
    colorArray18, &
    colorArray19 /)
#endif
end module m_wkml_color_def
