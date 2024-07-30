#compile helper functions
echo ========================================Building helpers=============================================
#EPI3V
/usr/bin/c++  -g -pg -O3 -w \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-I/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov \
-I/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK \
-I /home/jstewart23/Code/C++/Kapila \
-c  /home/jstewart23/Code/C++/Kapila/Epi3V.cpp -o Epi3V.o

echo =========================================Helpers built===============================================

echo ====================================Compiling Kapila.cpp================================
#export TP=~/Code/C++/TCHEM3.0Serial
#================
#Build ZeroDVar.o
#================
#Tons of depreciation warnings, suppressed
/usr/bin/c++ -w -O3  \
-I /home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov \
-I /home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK/ \
-I /home/jstewart23/Code/C++/Sundials_6_2/include/nvector/ \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-I /home/jstewart23/Code/C++/Kapila \
-c /home/jstewart23/Code/C++/Kapila/Kapila.cpp -o Kapila.o

#link in the new file: /home/jstewart23/Code/C++/Sundials3_0/cvode-6.1.1/installDir
#Void
echo ============================================Linking==================================================
#=======================================
#Link all the helpers into OneDIgnVar.o
#=======================================
#==============
#Modern version
#==============
/usr/bin/c++ Kapila.o \
-O3 -w \
~/Code/C++/Kapila/Epi3V.o \
/usr/lib/x86_64-linux-gnu/libdl.so \
/home/jstewart23/Code/C++/epic-cpp/build/Integrators/libepic1.0.0.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_generic.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_nvecserial.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_cvode.a \
/home/jstewart23/Code/C++/TChem/build/install/openblas/lib/libopenblas.so \
-o Kapila.x

echo =======================================Removing intermediary objects=================================

#Clean up
rm Kapila.o
./Kapila.x 1e-2 "CVODE" "Kapila.txt" 1e-10 1e-5 
./Kapila.x 1e-2 "EPI3V" "Kapila.txt" 1e-10 1e-1
