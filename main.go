package fft

import (
	"bytes"
	"math"
)

func log2(n int) int {
	return int(math.Log2(float64(n)))
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

func newFFTBuffer(n int) [][][]complex128 {
	k := log2(n)
	buff := make([][][]complex128, k+1)
	for i := k; i >= 0; i-- {
		buff[k-i] = make([][]complex128, n/(1<<i))
		for j := 0; j < n/(1<<i); j++ {
			buff[k-i][j] = make([]complex128, 1<<i)
		}
	}
	return buff
}

func nRootOfOne(n int) complex128 {
	return complex(math.Cos(2*math.Pi/float64(n)), math.Sin(2*math.Pi/float64(n)))
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

func waitChan(c *chan any, n int) {
	for i := 0; i < n; i++ {
		<-*c
	}
}

func fft(p *[][][]complex128, lvl int, v int, w complex128, c *chan any) {
	if len((*p)[lvl][v]) == 1 {
		*c <- nil
		return
	}
	n := len((*p)[lvl][v])
	k := n / 2
	for i := 0; i < n; i += 2 {
		(*p)[lvl+1][v*2][i/2] = (*p)[lvl][v][i]
		(*p)[lvl+1][v*2+1][i/2] = (*p)[lvl][v][i+1]
	}
	ch := make(chan any)
	go fft(p, lvl+1, v*2, w*w, &ch)
	go fft(p, lvl+1, v*2+1, w*w, &ch)
	waitChan(&ch, 2)
	var wt = complex(1, 0)
	for i := 0; i < n; i++ {
		(*p)[lvl][v][i] = (*p)[lvl+1][v*2][i%k] + (*p)[lvl+1][v*2+1][i%k]*wt
		wt *= w
	}
	*c <- nil
}

func evaluate(p []int) []complex128 {
	buff := newFFTBuffer(len(p))
	buff[0][0] = toComplex(p)
	ch := make(chan any)
	go fft(&buff, 0, 0, nRootOfOne(len(p)), &ch)
	waitChan(&ch, 1)
	return buff[0][0]
}

func interpolate(p []complex128) []int {
	buff := newFFTBuffer(len(p))
	buff[0][0] = p
	ch := make(chan any)
	go fft(&buff, 0, 0, nRootOfOne(-len(p)), &ch)
	waitChan(&ch, 1)
	res := toInt(buff[0][0])
	for i, val := range res {
		res[i] = val / len(p)
	}
	return res
}

type WideInt struct {
	val []int
	f   int
}

func (w *WideInt) size() int {
	return len(w.val)
}

func (w *WideInt) resize(n int) {
	for w.size() < n {
		w.val = append(w.val, 0)
	}
}

func (w *WideInt) carry() {
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

func (w *WideInt) Print() {
	if w.f == 1 {
		print("-")
	}
	var buf bytes.Buffer
	var i = w.size() - 1
	for ; i > 0 && w.val[i] == 0; i-- {
	}
	for ; i >= 0; i-- {
		buf.WriteByte(byte(w.val[i]) + '0')
	}
	print(buf.String())
}

func newWideInt(val []int, f int) WideInt {
	n := len(val)
	pow := 1
	for pow < n {
		pow <<= 1
	}
	resizedVar := make([]int, pow)
	copy(resizedVar, val)
	return WideInt{resizedVar, f}
}

func ToWideInt(val int) WideInt {
	var res = WideInt{[]int{abs(val)}, sign(val)}
	res.carry()
	return res
}

func ToInt(w WideInt) int {
	res := 0
	for i := w.size() - 1; i >= 0; i-- {
		res = res*10 + w.val[i]
	}
	return res
}

func Less(a WideInt, b WideInt) bool {
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

func Greater(a WideInt, b WideInt) bool {
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

func Equal(a WideInt, b WideInt) bool {
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

func LessOrEqual(a WideInt, b WideInt) bool {
	return Equal(a, b) || Less(a, b)
}

func GreaterOrEqual(a WideInt, b WideInt) bool {
	return Equal(a, b) || Greater(a, b)
}

func Add(a WideInt, b WideInt) WideInt {
	res := ToWideInt(0)
	if a.f == 1 && b.f == 1 {
		a.f, b.f = 0, 0
		res = Add(a, b)
		res.f = 1
		return res
	} else if a.f == 1 && b.f == 0 {
		a.f = 0
		return Subtract(b, a)
	} else if a.f == 0 && b.f == 1 {
		b.f = 0
		return Subtract(a, b)
	}
	sz := nextPowerOfTwo(max(a.size(), b.size()) + 1)
	a.resize(sz)
	b.resize(sz)
	res.resize(sz)
	for i := 0; i < sz; i++ {
		res.val[i] = (1-a.f*2)*a.val[i] + (1-b.f*2)*b.val[i]
	}
	for i := 0; i < sz-1; i++ {
		if res.val[i] < 0 {
			res.val[i] += 10
			res.val[i+1]--
		}
	}
	res.carry()
	return res
}

func Subtract(a WideInt, b WideInt) WideInt {
	res := ToWideInt(0)
	if b.f == 1 {
		b.f = 0
		return Add(a, b)
	} else if a.f == 1 {
		b.f = 1
		return Add(a, b)
	} else if a.f == 0 {
		if Less(a, b) {
			a.f, b.f = b.f, a.f
			res = Subtract(b, a)
			res.f = 1
			return res
		}
	}
	sz := nextPowerOfTwo(max(a.size(), b.size()) + 1)
	a.resize(sz)
	b.resize(sz)
	res.resize(sz)
	for i := 0; i < sz; i++ {
		res.val[i] = (1-a.f*2)*a.val[i] - (1-b.f*2)*b.val[i]
	}
	for i := 0; i < sz-1; i++ {
		if res.val[i] < 0 {
			res.val[i] += 10
			res.val[i+1]--
		}
	}
	res.carry()
	return res
}

func Multiply(a WideInt, b WideInt) WideInt {
	var sz = nextPowerOfTwo(a.size() + b.size() + 1)
	a.resize(sz)
	b.resize(sz)
	aDots := evaluate(a.val)
	bDots := evaluate(b.val)
	resDots := make([]complex128, sz)
	var res WideInt
	res.resize(sz)
	res.f = a.f ^ b.f
	for i := 0; i < sz; i++ {
		resDots[i] = aDots[i] * bDots[i]
	}
	res.val = interpolate(resDots)
	res.carry()
	return res
}

func Divide(a WideInt, b WideInt) WideInt {
	a.val = makeCopy(a.val)
	var base = 2
	var res = ToWideInt(0)
	var f = a.f ^ b.f
	a.f = 0
	b.f = 0
	for LessOrEqual(b, a) {
		var _b = newWideInt(b.val, b.f)
		var tmp = ToWideInt(1)
		for LessOrEqual(Multiply(_b, ToWideInt(base)), a) {
			_b = Multiply(_b, ToWideInt(base))
			tmp = Multiply(tmp, ToWideInt(base))
		}
		a = Subtract(a, _b)
		res = Add(res, tmp)
	}
	res.f = f
	return res
}

func Mod(a WideInt, b WideInt) WideInt {
	f := a.f
	a.f = 0
	return WideInt{Subtract(a, Multiply(b, Divide(a, b))).val, f}
}

func Pow(a WideInt, n WideInt) WideInt {
	var res = ToWideInt(1)
	for Greater(n, ToWideInt(0)) {
		if Equal(Mod(n, ToWideInt(2)), ToWideInt(1)) {
			res = Multiply(res, a)
		}
		a = Multiply(a, a)
		n = Divide(n, ToWideInt(2))
	}
	return res
}
