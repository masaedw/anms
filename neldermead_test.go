package neldermead

import (
	"math"
	"testing"
)

func norm(x Vec) float64 {
	n := 0.0
	for _, v := range x {
		n += math.Pow(v, 2)
	}
	return math.Sqrt(n)
}

const tol = 1e-8

func Test1(t *testing.T) {
	a := 1.0
	b := 100.0
	f := func(x Vec) float64 {
		return math.Pow(a-x[0], 2) + b*math.Pow(x[1]-math.Pow(x[0], 2), 2)
	}
	x0 := []float64{0, 0}
	x, fmin, err := Anms(f, x0)
	if err != nil {
		t.Error(err)
	}

	t1 := []float64{x[0] - a, x[1] - a}

	if norm(t1) >= tol {
		t.Errorf("norm too large %f", norm(t1))
	}

	if math.Abs(fmin-0.0) >= tol {
		t.Errorf("fmin too large %f", fmin)
	}
}
