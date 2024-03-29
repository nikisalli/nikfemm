# nikfemm
a 2d physics simulator written in C++ from scratch

# how to test it
clone the repo first:

``` git clone https://github.com/nikisalli/nikfemm.git ```

enter the repo

``` cd nikfemm ```

update submodules

``` git submodule update --init --remote ```

create build dir

``` mkdir build && cd build ```

configure and build

``` cmake .. ```

``` make -j4 ```

run an example

``` ../bin/test ```

# simulation example results

https://user-images.githubusercontent.com/31286021/214172922-65dab262-b9d2-44be-afaf-7bbbd28cfe62.mp4

![alt text](https://github.com/nikisalli/nikfemm/raw/main/images/iron.jpg "B plot iron C electromagnet and I iron piece")

![alt text](https://github.com/nikisalli/nikfemm/raw/main/images/halbach.jpg "halbach array")

![alt text](https://github.com/nikisalli/nikfemm/raw/main/images/motor.jpg "outrunner BLDC motor 2")

![alt text](https://github.com/nikisalli/nikfemm/raw/main/images/conductors.jpg "magnetic vector potential of rectangular conductors")

https://user-images.githubusercontent.com/31286021/214172986-e43e3fb3-43c8-4816-b682-b8a7dd3aaa1b.mp4
