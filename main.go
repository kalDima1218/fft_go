package main

import (
	"bytes"
	"fmt"
	"math"
)

var M = 998244353
var W = [24]int{1, 998244352, 911660635, 625715529, 373294451, 827987769, 280333251, 581015842, 628092333, 300892551, 586046298, 615001099, 318017948, 64341522, 106061068, 304605202, 631920086, 857779016, 841431251, 805775211, 390359979, 923521, 961, 31}

func makeCopy(a []int) []int {
	var res = make([]int, len(a))
	copy(res, a)
	return res
}

func sign(x int) int {
	if x < 0 {
		return 1
	}
	return 0
}

func log2(n int) int {
	return int(math.Log2(float64(n)))
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

func inv(a int) int {
	var n = M - 2
	var res = 1
	for n > 0 {
		if n&1 == 1 {
			res *= a
			res %= M
		}
		a *= a
		a %= M
		n /= 2
	}
	return res
}

func fft(p []int, w int) []int {
	if len(p) == 1 {
		return p
	}
	n := len(p)
	k := n / 2
	a := make([]int, k)
	b := make([]int, k)
	for i := 0; i < n; i += 2 {
		a[i/2] = p[i]
		b[i/2] = p[i+1]
	}
	a = fft(a, (w*w)%M)
	b = fft(b, (w*w)%M)
	wt := 1
	for i := 0; i < n; i++ {
		p[i] = (a[i%k] + b[i%k]*wt) % M
		wt = (wt * w) % M
	}
	return p
}

func evaluate(p []int) []int {
	return fft(p, W[log2(len(p))])
}

func interpolate(p []int) []int {
	p = fft(p, inv(W[log2(len(p))]))
	var invN = inv(len(p))
	for i, val := range p {
		p[i] = (val * invN) % M
	}
	return p
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

func swap(a *wideInt, b *wideInt) {
	*a, *b = *b, *a
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
	a.val = evaluate(makeCopy(a.val))
	b.val = evaluate(makeCopy(b.val))
	var res wideInt
	res.resize(sz)
	res.f = a.f ^ b.f
	for i := 0; i < sz; i++ {
		res.val[i] = (a.val[i] * b.val[i]) % M
	}
	res.val = interpolate(makeCopy(res.val))
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
