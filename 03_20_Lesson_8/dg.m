% approximation of the derivativve of g(gamma)

function y=dg(gamma)
epsilon=1e-7;
y=(g(gamma+epsilon)-g(gamma-epsilon))/(2*epsilon);