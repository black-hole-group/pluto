PLUTO v4.2
====================

![](http://astrologynewsservice.com/wp-content/uploads/2014/08/1240-305x260.png) 

Copyright (C) 2002-2015 Andrea Mignone. See [CONTRIBUTORS](./CONTRIBUTORS).

PLUTO is Godunov-type code for astrophysical fluid dynamics supporting several modules and algorithms. This is [version 4.2 (August 2015)](http://plutocode.ph.unito.it) of the code, with a few updates and bug fixes by the [Black Hole Group](https://blackholegroup.org).

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation. This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU General Public License](./LICENSE) for more details.

# System requirements

 - C compiler (gcc works fine),
 - Python (v. 2.0 or higher)

 Optionals

 - MPI (for parallel runs)
 - Chombo lib (for AMR)
 - libpng to produce on-the-fly graphics.
 

# Basic Installation and Usage

There's no `configure.sh`; once you've unpacked the distribution,  

1) define the shell variable `PLUTO_DIR` as the 
   main PLUTO directory, e.g., 
   if you are using tcsh:

```
setenv PLUTO_DIR "/home/user/PLUTO
```

   if you're using bash:

```
export PLUTO_DIR="/Users/ovidiu/PLUTO"
```

2) select a working dir anywhere on your hard disk; at the command prompt, just type 

```
python $PLUTO_DIR/setup.py
```

configure your problem and select makefile;

3) edit your `init.c` and `pluto.ini` to assign initial conditions and problem specific information;

4) compile with

```
make 
```

or `gmake`;

5) run with 

```
./pluto
```

See the documentation in `Doc/` for more information.
Have fun!
  
If you encounter problems during the previous steps, or have any other question regarding the code, please send an e-mail to `mignone@ph.unito.it`.



