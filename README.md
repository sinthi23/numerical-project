# Numerical Methods Laboratory

A comprehensive collection of numerical methods implementations in C++ for solving linear and non-linear equations.

---

# Table of Contents
- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)
- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton-Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)

---

## Solution of Linear Equations

A linear equation is an equation in which the highest power (degree) of the variable is 1.

**Example:**
- 3x - 2 = 0

---

### Gauss Elimination Method

#### Gauss Elimination Theory

[View Theory](./Linear%20Equation%20Methods/Gauss%20Elimination%20Method/Theory.md)

#### Gauss Elimination Code

[View Code](./Linear%20Equation%20Methods/Gauss%20Elimination%20Method/code.cpp)

#### Gauss Elimination Input

[View Input](./Linear%20Equation%20Methods/Gauss%20Elimination%20Method/input.txt)

#### Gauss Elimination Output

[View Output](./Linear%20Equation%20Methods/Gauss%20Elimination%20Method/output.txt)

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory

[View Theory](./Linear%20Equation%20Methods/Gauss%20Jordan%20Elimination%20Method/Theory.md)

#### Gauss Jordan Code

[View Code](./Linear%20Equation%20Methods/Gauss%20Jordan%20Elimination%20Method/code.cpp)

#### Gauss Jordan Input

[View Input](./Linear%20Equation%20Methods/Gauss%20Jordan%20Elimination%20Method/input.txt)

#### Gauss Jordan Output

[View Output](./Linear%20Equation%20Methods/Gauss%20Jordan%20Elimination%20Method/output.txt)

---

### LU Decomposition Method

#### LU Decomposition Theory

To factor any square matrix into two triangular matrices, i.e., one is a lower triangular matrix and the other is an upper triangular matrix, we can use the following steps.

**Steps for LU Decomposition:**

Start with a square matrix A: Given a square matrix A of size n ×n, the goal is to factor it into the product of two matrices: A = L×U, where:

- L is a lower triangular matrix with 1s on the diagonal.

- U is an upper triangular matrix.

Apply Gaussian elimination to convert matrix A into upper triangular form U. This step involves row operations to eliminate elements below the diagonal, resulting in an upper triangular matrix.

As we perform row operations, keep track of the multipliers used to eliminate the elements below the diagonal. These multipliers form the entries of the lower triangular matrix L.

- The entries of L will be the factors used during the elimination steps.

- The diagonal of L will consist of 1s.

#### LU Decomposition Code

[View Code](./Linear%20Equation%20Methods/LU%20Decomposition%20Method/code.cpp)

#### LU Decomposition Input

[View Input](./Linear%20Equation%20Methods/LU%20Decomposition%20Method/input.txt)

#### LU Decomposition Output

[View Output](./Linear%20Equation%20Methods/LU%20Decomposition%20Method/output.txt)

---

### Matrix Inversion

#### Matrix Inversion Theory

AX = B

A^(-1)AX = A^(-1)B

IX = A^(-1)B

X = A^(-1)B

**Steps to find the inverse of a matrix:**

1. Find the determinant.
2. Find the cofactor matrix.
3. Find the adjoint (Transpose of cofactor matrix).
4. Compute the inverse: A^(-1) = (1 / (|A|)) * adj(A).

#### Matrix Inversion Code

[View Code](./Linear%20Equation%20Methods/Matrix%20Inversion%20Method/code.cpp)

#### Matrix Inversion Input

[View Input](./Linear%20Equation%20Methods/Matrix%20Inversion%20Method/input.txt)

#### Matrix Inversion Output

[View Output](./Linear%20Equation%20Methods/Matrix%20Inversion%20Method/output.txt)

---

## Solution of Non-Linear Equations

A non-linear equation is one in which the power of the variable is greater than 1 or the variables appear in non-linear forms (squared, cubed, exponential, trigonometric, logarithmic, etc.).

**Examples:**
- e^x - 3x = 0
- x^2 - 4x - 10 = 0
- sin(x) = x/2

---

### Bisection Method

#### Bisection Theory

- Binary chopping or half-interval method

- If f(x) is real and continuous in the interval a < x < b, and f(a) and f(b) are of opposite sign, that is, f(a) f(b) < 0, then there is at least one real root in the interval between a and b.

- Let x1 = a and x2 = b. Define x0 to be the midpoint between a and b, that is,

X0 =(x1 + x2)/2

- Now there exist the following 3 conditions:
- (i) If f(x0) = 0, then the root is x0.
- (ii) If f(x0) * f(x1) < 0, then root is between x0 and x1.
- (iii) If f(x0) * f(x2) < 0, then root is between x0 and x2.

#### Bisection Code

[View Code](./Non-Linear%20Equation%20Methods/Bisection%20Method/code.cpp)

#### Bisection Input

[View Input](./Non-Linear%20Equation%20Methods/Bisection%20Method/input.txt)

#### Bisection Output

[View Output](./Non-Linear%20Equation%20Methods/Bisection%20Method/output.txt)

---

### False Position Method

#### False Position Theory

- Regular Falsi method in Latin

- Linear interpolation method

- If f(x) is real and continuous in the interval a < x < b, and f(a) and f(b) are of opposite sign, that is, f(a)f(b) < 0, then there is at least one real root in the interval between a and b.

- Let x1 = a and x2 = b. Define x0 to be the midpoint between a and b, that is,

X0 = x1 - f(x1) * (x2 – x1)/ (f(x2) - f(x1))

- Now there exist the following 3 conditions:
- (i) If f(x0) = 0, then the root is x0.
- (ii) If f(x0) * f(x1) < 0, then root is between x0 and x1.
- (iii) If f(x0) * f(x2) < 0, then root is between x0 and x2.

#### False Position Code

[View Code](./Non-Linear%20Equation%20Methods/False%20Position%20Method/code.cpp)

#### False Position Input

[View Input](./Non-Linear%20Equation%20Methods/False%20Position%20Method/input.txt)

#### False Position Output

[View Output](./Non-Linear%20Equation%20Methods/False%20Position%20Method/output.txt)

---

### Newton-Raphson Method

#### Newton-Raphson Theory

- If f(x) is a real and continuously differentiable function, the Newton–Raphson method is used to find a root of the equation f(x) = 0.

- Let the initial approximation to the root be x0.

- A better approximation x1 is obtained using the formula: x1 = x0 – (f(x0)/f'(x0)).

- This process is repeated iteratively as:

X(n+1) = xn – (f(xn)/f'(xn))

- The iteration continues until the difference between two successive values of x becomes sufficiently small.

#### Newton-Raphson Code

[View Code](./Non-Linear%20Equation%20Methods/Newton-Raphson%20Method/code.cpp)

#### Newton-Raphson Input

[View Input](./Non-Linear%20Equation%20Methods/Newton-Raphson%20Method/input.txt)

#### Newton-Raphson Output

[View Output](./Non-Linear%20Equation%20Methods/Newton-Raphson%20Method/output.txt)

---

### Secant Method

#### Secant Theory

Secant method is a recursive method for finding the root of a polynomial by successive approximation. In this method, the neighborhood roots are approximated by secant line (A line that intersects the curve at two distinct points) to the function f(x).

Formula for Secant Method: xk + 1 = xk – (((xk – 1) – xk)/ (f(xk-1) – f(xk))*f(xk)

The secant method will find a root of a function f if the starting points x0 and x1 are close enough to the actual root and if the function behaves well. If the function is smooth and has a simple root (i.e., the root only occurs once), the method converges to the root at a rate related to the golden ratio: Φ = 1.618

This means the method improves quickly but not as fast as some other methods.

#### Secant Code

[View Code](./Non-Linear%20Equation%20Methods/Secant%20Method/code.cpp)

#### Secant Input

[View Input](./Non-Linear%20Equation%20Methods/Secant%20Method/input.txt)

#### Secant Output

[View Output](./Non-Linear%20Equation%20Methods/Secant%20Method/output.txt)

---

## Repository Structure
