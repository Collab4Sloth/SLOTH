# Multilinear interpolation

This test validates a class for multivariate interpolation in dimension n.

Trilinear interpolation can be written as :

$$	
f(x,y,z) = (1 - t)  (1 - u) (1-v) f_{i,j,k} + t (1-u) (1-v) f_{i+1,j,k} + u (1-t) (1-v) f_{i,j+1,k} + v (1-t) (1-u) f_{i,j,k+1} \\
	 + t u (1-v) f_{i+1,j+1,k} + t v (1-u) f_{i+1,j,k+1} + u v (1-t)f_{i,j+1,k+1} + t u v f_{i+1,j+1,k+1} 
$$
this formula can be generalized as :
$$
f(\bm{x}) = \sum_{i=0}^{2^N-1}\left(\prod_{k=1}^{N}w_{i,k}\right)f(\bm{\tilde{x}})
$$
where
\[
w_{i,k} = \left\{
\begin{array}{ll}
1 - \alpha_i    &\text{if } \bm{\tilde{x}} \text{ uses the downstream point} \\
\alpha_i        &\text{if } \bm{\tilde{x}} \text{ uses the upstream point}
\end{array}
\right.
\]
and 
$$
\alpha_i = \frac{x^i - \tilde{x}^i_j}{\tilde{x}^i_{j+1}- \tilde{x}^i_j}
$$
  

## Running

### __Using the binary__
```shell
./UtilitiesMultilinearInterpolator
```

### __Using ctest__

```shell
ctest -R UtilitiesMultilinearInterpolator
```

### __In case of code coverage analysis__

```shell
make UtilitiesMultilinearInterpolator
```


## Post-processing

None

## Files & Dependencies

Source file : `main.cpp`

## References

None

## Intellectual Property

See [About page](../../../../../about.html) 


