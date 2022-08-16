
%chk=gaussian
%rwf=gaussian
%nproc=6
%mem=8GB

#P SP PBE1PBE/6-311G* SCF=VeryTight Integral(Grid=UltraFine) Force NoSymm

DFT single-point calculation

0 1
C                -0.610677   -0.358946    0.000000
O                -1.722641    0.103619    0.000000
C                 0.622629    0.441868    0.000000
C                 1.819262   -0.142900    0.000000
H                -0.448906   -1.460689    0.000000
H                 0.501288    1.521894    0.000000
H                 1.913409   -1.226676    0.000000
H                 2.743326    0.425430    0.000000



