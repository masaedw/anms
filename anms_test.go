package anms

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

func Test2(t *testing.T) {
	// asset allocation re-balancing without selling

	current := []float64{300, 200, 400, 100}
	targetRatio := []float64{0.25, 0.25, 0.25, 0.25}

	n := 500.0
	total := 300.0 + 200 + 400 + 100 + n
	targetNum := make([]float64, len(current))
	for i := range targetNum {
		targetNum[i] = total * targetRatio[i]
		t.Logf("targetNum[%d] = %f\n", i, targetNum[i])
	}

	f := func(x Vec) float64 {
		r := 0.0
		for i := range x {
			r += math.Pow(targetNum[i]-current[i]-x[i], 2)
		}

		sum := 0.0
		for _, v := range x {
			sum += v
			if v < 0 {
				return math.Inf(1)
			}
		}

		if n < sum {
			return math.Inf(1)
		}

		return r + math.Pow(sum-n, 16)
	}

	x0 := []float64{0, 0, 100, 0}
	x, fmin, err := Anms(f, x0)
	if err != nil {
		t.Error(err)
	}

	t.Log(x)
	t.Log(fmin)
	t.Log(x[0] + x[1] + x[2] + x[3])
}
