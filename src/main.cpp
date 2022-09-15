#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv){
    cout << "You have entered " << argc
         << " arguments:" << "\n";
  
    for (int i = 0; i < argc; ++i)
        cout << argv[i] << "\n";
  
    return 0;
}