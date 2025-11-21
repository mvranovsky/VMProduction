# Vector Meson production

Code used for calculation of vector meson production in UPC collisions or diffractive proton-proton collisions  .

Before compiling, check you have installed all the dependencies: GSL, Cuba, and ROOT for plotting.

To compile, use the commands below:

```bash
mkdir build
cd build 
```

which will create a *build/* directory. Inside the directory, use the following commands to compile 

```bash
cmake ../
make
```

Then, just add source codes to *src/* directory and header files to *include/*. Whenever compiling again, use only 

```bash
make
```

After successful compilation, file *build/VMProduction* executable will appear. Run it from *build/* dir

```bash
./VMProduction
```
