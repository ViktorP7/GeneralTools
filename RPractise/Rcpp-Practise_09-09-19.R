# Practice writing C++ in R using https://adv-r.hadley.nz/rcpp.html

# Call package
library(Rcpp, bench)

# Create simple addition function
cppFunction('int add(int x, int y, int z){
  int sum = x + y + z;
  return sum;
}')

add(2,5,7)

# Create the sign function in cpp, which returns 1 if input positive, -1 if negative
cppFunction('int signC(int x) {
              if (x > 0) {
                return 1;
              } else if (x == 0) {
                return 0;
              } else {
                return -1;
              }
            }')

signC(-664)
            
# Create sum function in cpp - loops in cpp are much more efficient than R
cppFunction('double sumC(NumericVector x) {
              /* This calculates the size of the input vector x */
              int n = x.size(); 
              /* Initialise the total variable */
              double total = 0;
              /* C++ for loop has syntax for(initialise; check; increment)
              * Initialise variable i=0, check if i<n (loop stops if not), increment by 1 using ++
              * C++ indices start at 0, so i=0 */
              for(int i = 0; i < n; ++i) {
                /* C++ has special operators such as += that modify in place */
                total += x[i];
              }
              return total;
            }')
v <- c(1,2,3,4,5)

sumC(v)

# Create a function to calculate euclidean distance between a point and vector, returning a vector
cppFunction('NumericVector pdistC(double x, NumericVector ys) {
            
              int n = ys.size();

              /* Create numeric vector of length n */            
              NumericVector out(n);
            
              /* cpp uses pow() instead of ^ for exponents */
              for(int i = 0; i < n; ++i) {
                out[i] = sqrt(pow(ys[i] - x, 2.0));
              }
            
              return out;
            }')

pdistC(4, v)

# Call the cpp script with the needed function
pathCpp <- "C://Users/UCD/Desktop/UbuntuSharedFolder/GeneralTools/CppPractise/Cpp-MeanFunction_09-09-19.cpp"
sourceCpp(pathCpp)

meanC(v)

# Call cpp script with task functions
