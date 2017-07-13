# UODESys

An open-source MATLAB toolbox for the reduction of high-dimensional ODE systems to a low-dimensional, uncertain system.

## Contents
- [System requirements](#Requirements)
- [License](#License)
- [How to use](#HowToUse)

## System requirements<a name="Requirements"></a>

In order to use QUINOPT, you will need:

1. A working version of [YALMIP](https://yalmip.github.io/), the MATLAB optimization modelling software by J. L&ouml;fberg
2. A suitable SDP solver. Choices include [SeDuMi](https://github.com/sqlp/sedumi), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [SDPA](http://sdpa.sourceforge.net/), [Mosek](https://www.mosek.com/) (free for
    users in academia).

## License<a name="License"></a>

UODESys is distributed under the [GNU GENERAL PUBLIC LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html).

Copyright 2017, M. Lakshmi.

Licensed under the GNU General Public License, Version 3.0 (the "License"); you may not use this file except in compliance with the License. You should have received a copy of the License with this repository, in a file called LICENSE.txt. If this is not the case, you may obtain a copy of the License at [https://www.gnu.org/licenses/gpl-3.0.en.html](https://www.gnu.org/licenses/gpl-3.0.en.html).

**Disclaimer of warranty:** There is no warranty for the program, to the extent permitted by
applicable law.  Except when otherwise stated in writing the copyright
holders and/or other parties provide the program "as is" without warranty
of any kind, either expressed or implied, including, but not limited to,
the implied warranties of merchantability and fitness for a particular
purpose.  The entire risk as to the quality and performance of the program
is with you.  Should the program prove defective, you assume the cost of
all necessary servicing, repair or correction.

**Limitation of liability:** In no event unless required by applicable law or agreed to in writing
will any copyright holder, or any other party who modifies and/or conveys
the program as permitted above, be liable to you for damages, including any
general, special, incidental or consequential damages arising out of the
use or inability to use the program (including but not limited to loss of
data or data being rendered inaccurate or losses sustained by you or third
parties or a failure of the program to operate with any other programs),
even if such holder or other party has been advised of the possibility of
such damages.

## How to use UODESys<a name="HowToUse"></a>

To install UODESys, simply download or clone this repository and add all files to your MATLAB path. For a quick introduction to UODESys, please refer to the user manual, which you can find inside the `docs/` folder.
