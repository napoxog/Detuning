#!/usr/bin/python
__author__ = 'agritsenko'
import os
import sys

import numpy as np


# spline interpolation procedure
def spl_int(x, y, x_out, y_out, step):
    n = len(x)
    x0 = x[0]
    xN = x[n - 1]
    b = x.copy()
    c = x.copy()
    d = x.copy()
    old_step = (xN - x0) / n
    # print ("len=",n,"range=",x0,xN,step)
    x_out = np.arange(x0, xN, step)
    y_out = x_out.copy()
    # print ("new_len=",len(x_out),"old_step=",old_step)

    #   calc spline koef
    if n == 2:
        b[0] = b[1] = (y[1] - y[0]) / (x[1] - x[0])
        c[0] = c[1] = d[0] = d[1] = 0
    else:
        d[0] = x[1] - x[0]
        c[1] = (y[1] - y[0]) / d[0]
        for i in range(1, n - 1):
            d[i] = x[i + 1] - x[i]
            b[i] = 2.0 * (d[i - 1] + d[i])
            c[i + 1] = (y[i + 1] - y[i]) / d[i]
            c[i] = c[i + 1] - c[i]
        b[0] = -d[0]
        b[n - 1] = -d[n - 2]
        c[0] = c[n - 1] = 0
        if n != 3:
            c[0] = (c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0])) * d[0] * d[0] / (x[3] - x[0]);
            c[n - 1] = -(c[n - 2] / (x[n - 1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4])) * d[n - 2] * d[n - 2] / (
                    x[n - 1] - x[n - 4]);
        for i in range(1, n):
            b[i] = b[i] - d[i - 1] * d[i - 1] / b[i - 1]
            c[i] = c[i] - d[i - 1] * c[i - 1] / b[i - 1]
        c[n - 1] /= b[n - 1]
        for j in range(1, n):
            i = n - j - 1
            c[i] = (c[i] - d[i] * c[i + 1]) / b[i]
        b[n - 1] = (y[n - 1] - y[n - 2]) / d[n - 2] + d[n - 2] * c[n - 2] + 2 * c[n - 1]
        for i in range(0, n - 1):
            b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2. * c[i])
            d[i] = (c[i + 1] - c[i]) / d[i]
            c[i] = 3.0 * c[i]
        c[n - 1] = 3.0 * c[n - 1]
        d[n - 1] = d[n - 2]
    # /*    calc interpolation splain */

    lx = x[0]
    rx = x[n - 1]
    i = 0
    nmax = len(x_out)

    for j in range(0, nmax):
        # ; lx <= x[n-1], num < max(num_out_points,n) ; num=++j,lx+=new_step )
        r = x_out[j] - x[i]
        k = i + 1
        y_out[j] = y[i] + r * (b[i] + r * (c[i] + r * d[i]))
        # print (i,x[i],y[i],j,x_out[j],y_out[j],r)
        if (x[i + 1] < x_out[j]):
            i = i + 1

    #    y_out[j]= x[n-1]

    # return(outarr);
    return y_out.copy(), x_out.copy()


class Wavelet:
    name = ''
    file = ''
    Length = float(0)
    Samples = int(0)
    Samples_interp = int(0)
    SampleRate = float(0)
    SampleRate_interp = float(0.5)
    min = float(0)
    max = float(0)
    time = np.array(0.0)
    amp = np.array(0.0)
    time_interp = np.array(0.0)
    amp_interp = np.array(0.0)
    zeroIdx = int(0)
    zeroIdx_interp = int(0)
    min_idx = int(0)
    min_idx_interp = int(0)

    def Load(self, path):
        self.file = path
        data = [geom_line.strip() for geom_line in open(path)]
        line_idx = -1
        for line in data:
            if line_idx < 0:
                self.zeroIdx = float(line.split()[1])
                self.SampleRate = float(line.split()[2])
                self.Samples = float(line.split()[0])
                self.Length = self.Samples * self.SampleRate
            else:
                amp = float(line.split()[0])
                self.time = np.append(self.time, self.SampleRate * line_idx - self.zeroIdx * self.SampleRate)
                self.amp = np.append(self.amp, amp)

                # print self.amp.__len__()
                if self.min > amp:
                    self.min = amp
                    self.min_idx = abs(line_idx - w.zeroIdx)
                if self.max < amp:
                    self.max = amp
            # print (line_idx,line_idx*self.SampleRate-self.zeroIdx*self.SampleRate,line)
            line_idx = line_idx + 1
        self.time = np.delete(self.time, 0)
        self.amp = np.delete(self.amp, 0)


ScaleFactor = 1

if len(sys.argv) < 2:
    print("Please provide Wavelet ASCII File\n\n   Usage: Detuning.py Wavelet_file [ScaleFactor] [resampling_rate])")
    sys.exit(1)

if os.path.isfile(sys.argv[1]):
    w = Wavelet()
    w.Load(sys.argv[1])
else:
    print("Cannot find file:", sys.argv[1])
    sys.exit(2)

if len(sys.argv) > 2:
    ScaleFactor = float(sys.argv[2])

if len(sys.argv) > 3:
    w.SampleRate_interp = float(sys.argv[3])

# print ('length=',format(w.amp.__len__()))
print('scale_factor=', format(ScaleFactor))
print('sample_rate=', format(w.SampleRate))
print('resample_rate=', format(w.SampleRate_interp))
print('min_idx=', format(w.min_idx))
print('zero_idx=', format(w.zeroIdx))
print('min=', format(w.min))
print('max=', format(w.max))

# perform spline interpolation.
spline = spl_int(w.time, w.amp, w.time_interp, w.amp_interp, w.SampleRate_interp)
w.amp_interp = spline[0].copy()
w.time_interp = spline[1].copy()

# print w.time_interp[0],w.amp_interp[0]

w.zeroIdx_interp = w.zeroIdx * (float(w.SampleRate) / w.SampleRate_interp)
w.min_idx_interp = w.min_idx * (w.SampleRate / w.SampleRate_interp)
w.Samples_interp = len(w.amp_interp)

# exit()
print('min_idx_interp=', format(w.min_idx_interp))
print('zero_idx_interp=', format(w.zeroIdx_interp))

tuning_samples = w.Samples_interp / 2
tuning_length = tuning_samples * w.SampleRate_interp

tuning_values = []
tuning_segments = []

print('Samples=', tuning_samples, "length=", tuning_length)

# calcullating tuning curve and position of its max.
dt = tuning_length
dtidx = int(tuning_samples)
camp_max_idx = int(0)
camp_max = 0
while dtidx >= 0:
    zero_amp = w.amp_interp[int(w.zeroIdx_interp)]
    widx = int(max(w.zeroIdx_interp - dtidx, 0))
    a1 = zero_amp - w.amp_interp[widx]

    # print ("dt=",dt,"zero_amp=",zero_amp,"widx1=",widx)
    widx = int(min(w.zeroIdx_interp + dtidx, w.Samples_interp - 1))
    a2 = -zero_amp + w.amp_interp[widx]

    # print ("dt=",dt,"zero_amp=",zero_amp,"widx2=",widx)

    camp = ScaleFactor * (abs(a1) + abs(a2))
    # print dt,camp,a1,a2
    if camp > camp_max:
        camp_max = camp
        camp_max_idx = dtidx
    tuning_values.insert(0, (dt, camp, 0, 0))

    dt = dt - w.SampleRate_interp
    dtidx = dtidx - 1
# exit()
print('camp_max_idx=', format(camp_max_idx))

# calcullating no tune curve and transfer function
dt = tuning_length
dtidx = int(tuning_samples)

print('{0:20} {1:20}  {2:20} {3:20}'.format("delta_t", "Tuning_Curve", "Normal_Curve", "Transfer_Func"))
while dtidx >= 0:

    camp = tuning_values[dtidx][1]

    # print (tuning_values.__len__())
    notune = (tuning_values[len(tuning_values) - 1][1])
    if dtidx < camp_max_idx:
        # notune=math.atan(notune/float(camp_max_idx)/ScaleFactor)*dtidx*ScaleFactor
        notune = notune / camp_max * camp
        # print (format(notune/camp_max_idx),format(notune),format(camp_max_idx))

    if tuning_values[dtidx][0] != 0:
        tuning_values[dtidx] = (
            tuning_values[dtidx][0], tuning_values[dtidx][1], notune, notune / tuning_values[dtidx][0])
    else:
        tuning_values[dtidx] = (tuning_values[dtidx][0], tuning_values[dtidx][1], notune, 0)

    if dtidx:
        print('{0:10}  {1:20} {2:20} {3:20} {4:20} '.format(dt, camp, notune, notune / camp,
                                                            (tuning_values[len(tuning_values) - 1][1]) * notune / camp))
    else:
        print('{0:10}  {1:20} {2:20} {3:20} {4:20} '.format(dt, camp, notune, prev_notune / prev_camp,
                                                            prev_notune / prev_camp *
                                                            tuning_values[len(tuning_values) - 1][1]))
    dt = dt - w.SampleRate_interp
    prev_camp = camp
    prev_notune = notune
    dtidx = dtidx - 1
