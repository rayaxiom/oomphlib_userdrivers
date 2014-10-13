#include <cassert>

/// Start of fract0
// A recursive implementation of factorial
int factorial(int i)
{
    if(i)
        return i * factorial(i-1);
    else
        return 1;
}
/// End of fract0

//----------------main0----------------
// Test the implementation.
int main()
{
    assert(factorial(5) == 120);
}
//----------------main1----------------

