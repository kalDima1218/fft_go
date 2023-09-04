package main

import (
	"bytes"
	"fmt"
	"math"
)

func makeCopy(a []int) []int {
	var res = make([]int, len(a))
	copy(res, a)
	return res
}

func toComplex(a []int) []complex128 {
	var res = make([]complex128, len(a))
	for i := 0; i < len(a); i++ {
		res[i] = complex(float64(a[i]), 0)
	}
	return res
}

func toInt(a []complex128) []int {
	var res = make([]int, len(a))
	for i := 0; i < len(a); i++ {
		res[i] = int(math.Ceil(real(a[i])))
	}
	return res
}

func nRootOfOne(n int) complex128 {
	return complex(math.Cos(2*math.Pi/float64(n)), math.Sin(2*math.Pi/float64(n)))
}

func sign(x int) int {
	if x < 0 {
		return 1
	}
	return 0
}

func max(a int, b int) int {
	if a > b {
		return a
	}
	return b
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func nextPowerOfTwo(n int) int {
	n--
	n |= n >> 1
	n |= n >> 2
	n |= n >> 4
	n |= n >> 8
	n |= n >> 16
	n++
	return n
}

func fft(p []complex128, w complex128) []complex128 {
	if len(p) == 1 {
		return p
	}
	n := len(p)
	k := n / 2
	a := make([]complex128, k)
	b := make([]complex128, k)
	for i := 0; i < n; i += 2 {
		a[i/2] = p[i]
		b[i/2] = p[i+1]
	}
	a = fft(a, w*w)
	b = fft(b, w*w)
	var wt = complex(1, 0)
	for i := 0; i < n; i++ {
		p[i] = a[i%k] + b[i%k]*wt
		wt *= w
	}
	return p
}

func evaluate(p []int) []complex128 {
	return fft(toComplex(p), nRootOfOne(len(p)))
}

func interpolate(p []complex128) []int {
	res := toInt(fft(p, nRootOfOne(-len(p))))
	for i, val := range res {
		res[i] = val / len(p)
	}
	return res
}

type wideInt struct {
	val []int
	f   int
}

func (w *wideInt) size() int {
	return len(w.val)
}

func (w *wideInt) resize(n int) {
	for w.size() < n {
		w.val = append(w.val, 0)
	}
}

func (w *wideInt) carry() {
	size := w.size()
	for i := 0; i < size-1; i++ {
		w.val[i+1] += w.val[i] / 10
		w.val[i] %= 10
	}
	for abs(w.val[size-1]) > 9 {
		w.val = append(w.val, w.val[size-1]/10)
		size++
		w.val[size-2] %= 10
	}
	for size > 1 && w.val[size-1] == 0 {
		w.val = w.val[:size-1]
		size--
	}
	w.resize(nextPowerOfTwo(size))
}

func (w *wideInt) print() {
	if w.f == 1 {
		fmt.Print("-")
	}
	var buf bytes.Buffer
	var i = w.size() - 1
	for ; i > 0 && w.val[i] == 0; i-- {
	}
	for ; i >= 0; i-- {
		buf.WriteByte(byte(w.val[i]) + '0')
	}
	fmt.Println(buf.String())
}

func newWideInt(val []int, f int) wideInt {
	n := len(val)
	pow := 1
	for pow < n {
		pow <<= 1
	}
	resizedVar := make([]int, pow)
	copy(resizedVar, val)
	return wideInt{resizedVar, f}
}

func toWideInt(val int) wideInt {
	var res = wideInt{[]int{abs(val)}, sign(val)}
	res.carry()
	return res
}

func less(a wideInt, b wideInt) bool {
	if a.f != b.f {
		return a.f == 1
	}
	var sz = nextPowerOfTwo(max(a.size(), b.size()))
	a.resize(sz)
	b.resize(sz)
	for i := sz - 1; i >= 0; i-- {
		if a.val[i] < b.val[i] {
			return a.f != 1
		}
		if a.val[i] > b.val[i] {
			return a.f == 1
		}
	}
	return a.f == 1
}

func greater(a wideInt, b wideInt) bool {
	if a.f != b.f {
		return a.f != 1
	}
	var sz = nextPowerOfTwo(max(a.size(), b.size()))
	a.resize(sz)
	b.resize(sz)
	for i := sz - 1; i >= 0; i-- {
		if a.val[i] > b.val[i] {
			return a.f != 1
		}
		if a.val[i] < b.val[i] {
			return a.f == 1
		}
	}
	return a.f == 1
}

func equal(a wideInt, b wideInt) bool {
	if a.f != b.f {
		return false
	}
	var sz = nextPowerOfTwo(max(a.size(), b.size()))
	a.resize(sz)
	b.resize(sz)
	for i := sz - 1; i >= 0; i-- {
		if a.val[i] != b.val[i] {
			return false
		}
	}
	return true
}

func add(a wideInt, b wideInt) wideInt {
	sz := nextPowerOfTwo(max(a.size(), b.size()) + 1)
	a.resize(sz)
	b.resize(sz)
	res := toWideInt(0)
	res.resize(sz)
	for i := 0; i < sz; i++ {
		res.val[i] = (1-a.f*2)*a.val[i] + (1-b.f*2)*b.val[i]
	}
	res.carry()
	sz = res.size()
	for i := 0; i < sz-1; i++ {
		if res.val[i] < 0 {
			res.val[i] += 10
			res.val[i+1]--
		}
	}
	for i := 0; i < sz; i++ {
		if res.val[i] < 0 {
			res.val[i] = -res.val[i]
			res.f = (res.f + 1) % 2
			break
		}
	}
	return res
}

func subtract(a wideInt, b wideInt) wideInt {
	sz := nextPowerOfTwo(max(a.size(), b.size()) + 1)
	a.resize(sz)
	b.resize(sz)
	res := toWideInt(0)
	res.resize(sz)
	for i := 0; i < sz; i++ {
		res.val[i] = (1-a.f*2)*a.val[i] - (1-b.f*2)*b.val[i]
	}
	res.carry()
	sz = res.size()
	for i := 0; i < sz-1; i++ {
		if res.val[i] < 0 {
			res.val[i] += 10
			res.val[i+1]--
		}
	}
	for i := 0; i < sz; i++ {
		if res.val[i] < 0 {
			res.val[i] = -res.val[i]
			res.f = (res.f + 1) % 2
			break
		}
	}
	return res
}

func multiply(a wideInt, b wideInt) wideInt {
	var sz = nextPowerOfTwo(a.size() + b.size() + 1)
	a.resize(sz)
	b.resize(sz)
	aDots := evaluate(a.val)
	bDots := evaluate(b.val)
	resDots := make([]complex128, sz)
	var res wideInt
	res.resize(sz)
	res.f = a.f ^ b.f
	for i := 0; i < sz; i++ {
		resDots[i] = aDots[i] * bDots[i]
	}
	res.val = interpolate(resDots)
	res.carry()
	return res
}

func divide(a wideInt, b wideInt) wideInt {
	a.val = makeCopy(a.val)
	var base = 2
	var res = toWideInt(0)
	var f = a.f ^ b.f
	a.f = 0
	b.f = 0
	for less(b, a) || equal(b, a) {
		var _b = newWideInt(b.val, b.f)
		var tmp = toWideInt(1)
		for less(multiply(_b, toWideInt(base)), a) || equal(multiply(_b, toWideInt(base)), a) {
			_b = multiply(_b, toWideInt(base))
			tmp = multiply(tmp, toWideInt(base))
		}
		a = subtract(a, _b)
		res = add(res, tmp)
	}
	res.f = f
	return res
}

func mod(a wideInt, b wideInt) wideInt {
	return subtract(a, multiply(b, divide(a, b)))
}

func pow(a wideInt, n wideInt) wideInt {
	var res = toWideInt(1)
	for greater(n, toWideInt(0)) {
		if equal(mod(n, toWideInt(2)), toWideInt(1)) {
			res = multiply(res, a)
		}
		a = multiply(a, a)
		n = divide(n, toWideInt(2))
	}
	return res
}

func isPrime(n wideInt) bool {
	if less(n, toWideInt(2)) {
		return false
	}
	if equal(n, toWideInt(2)) {
		return true
	}
	if equal(mod(n, toWideInt(2)), toWideInt(0)) {
		return false
	}
	var i = toWideInt(3)
	for less(multiply(i, i), n) || equal(multiply(i, i), n) {
		if equal(mod(n, i), toWideInt(0)) {
			return false
		}
		i = add(i, toWideInt(2))
	}
	return true
}

func main() {
	var i = toWideInt(1)
	for true {
		if isPrime(i) {
			i.print()
		}
		i = add(i, toWideInt(1))
	}
}
