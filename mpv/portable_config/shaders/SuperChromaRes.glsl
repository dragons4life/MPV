// SuperChromaRes by Shiandow
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library.


//!HOOK NATIVE
//!BIND HOOKED
//!BIND CHROMA
//!SAVE DOWNSCALEDCHX
//!WIDTH CHROMA.w
//!COMPONENTS 4

// -- Downscaling --
#define offset   (-vec2(0.0, 0.0)*NATIVE_size*CHROMA_pt) // -chromaOffset*LUMA_size/CHROMA_size

#define dxdy     (vec2(CHROMA_pt.x, NATIVE_pt.y)) // 1/output_size
#define ddxddy   (NATIVE_pt)                      // 1/input_size

#define factor   ((ddxddy*vec2(CHROMA_size.x, NATIVE_size.y))[axis])

#define axis 0

#define Kernel(x) clamp(0.5 + (0.5 - abs(x)) / factor, 0.0, 1.0)
#define taps (1.0 + factor)

vec4 hook() {
    //if (CHROMA_size.x >= NATIVE_size.x) return vec4(0);
    // Calculate bounds
    float low  = floor((NATIVE_pos - 0.5*taps*dxdy) * NATIVE_size - offset + 0.5)[axis];
    float high = floor((NATIVE_pos + 0.5*taps*dxdy) * NATIVE_size - offset + 0.5)[axis];

    float W = 0.0;
    vec2 avg = vec2(0);
    vec2 pos = NATIVE_pos;

    for (int k = 0; k < int(high - low); k++) {
        pos[axis] = ddxddy[axis] * (float(k) + low + 0.5);
        float rel = (pos[axis] - NATIVE_pos[axis])*vec2(CHROMA_size.x, NATIVE_size.y)[axis] + offset[axis]*factor;
        float w = Kernel(rel);

        avg += w * textureLod(NATIVE_raw, pos, 0.0).yz;
        W += w;
    }
    avg /= vec2(W);

    return vec4(0.0, avg, 1.0);
}

//!HOOK NATIVE
//!BIND CHROMA
//!BIND DOWNSCALEDCHX
//!SAVE LOWRES_UV
//!WIDTH CHROMA.w
//!HEIGHT CHROMA.h
//!COMPONENTS 4

// -- Downscaling --
#define offset   (-vec2(0.0, 0.0)*DOWNSCALEDCHX_size*CHROMA_pt)

#define dxdy     (CHROMA_pt)
#define ddxddy   (DOWNSCALEDCHX_pt)

#define factor   ((ddxddy*CHROMA_size)[axis])

#define axis 1

#define Kernel(x) clamp(0.5 + (0.5 - abs(x)) / factor, 0.0, 1.0)
#define taps (1.0 + factor)

vec4 hook() {
    // Calculate bounds
    float low  = floor((DOWNSCALEDCHX_pos - 0.5*taps*dxdy) * DOWNSCALEDCHX_size - offset + 0.5)[axis];
    float high = floor((DOWNSCALEDCHX_pos + 0.5*taps*dxdy) * DOWNSCALEDCHX_size - offset + 0.5)[axis];

    float W = 0.0;
    vec4 avg = vec4(0);
    vec2 pos = DOWNSCALEDCHX_pos;

    for (int k = 0; k < int(high - low); k++) {
        pos[axis] = ddxddy[axis] * (float(k) + low + 0.5);
        float rel = (pos[axis] - DOWNSCALEDCHX_pos[axis])*CHROMA_size[axis] + offset[axis]*factor;
        float w = Kernel(rel);

        avg += w * textureLod(DOWNSCALEDCHX_raw, pos, 0.0);
        W += w;
    }
    avg /= vec4(W);

    return vec4(0.0, avg.yz, avg.x);
}

//!HOOK NATIVE
//!BIND HOOKED
//!BIND CHROMA
//!BIND LOWRES_UV

//	-- SuperChromaRes --

#define strength  1.0

#define acuity 100.0
#define radius 0.5
#define power 0.5

#define dxdy (HOOKED_pt)
#define ddxddy (LOWRES_UV_pt)
#define chromaOffset vec2(0.0, 0.0)

// -- Window Size --
#define taps 3.0
#define even (taps - 2.0 * floor(taps / 2.0) == 0.0)
#define minX int(1.0-ceil(taps/2.0))
#define maxX int(floor(taps/2.0))

#define sqr(x) dot(x,x)
#define factor (ddxddy*HOOKED_size)
#define Kernel(x) (cos(acos(-1.0)*(x)/taps)) // Hann kernel

//Current high res value
#define Get(x,y)    ( textureLod(HOOKED_raw, HOOKED_pos + dxdy*vec2(x,y), 0.0).xyz )
#define GetW(x,y)   ( textureLod(LOWRES_UV_raw, ddxddy*(pos+vec2(x,y)+0.5), 0.0).w )
//Downsampled result
#define LowRes(x,y) ( textureLod(LOWRES_UV_raw, ddxddy*(pos+vec2(x,y)+0.5), 0.0) )

vec4 hook() {    
    vec4 c0 = HOOKED_tex(HOOKED_pos);
    //if (LOWRES_UV_size.x >= HOOKED_size.x) return c0;

    // Calculate position
    vec2 pos = HOOKED_pos * LOWRES_UV_size - chromaOffset - vec2(0.5);
    vec2 offset = pos - (even ? floor(pos) : round(pos));
    pos -= offset;

    // Calculate faithfulness force
    float weightSum = 0.0;
    vec4 diff = vec4(0);

    for (int X = minX; X <= maxX; X++)
    for (int Y = minX; Y <= maxX; Y++)
    {
        float dI2 = sqr(acuity*(c0.x - GetW(X,Y)));
        //float dXY2 = sqr((vec2(X,Y) - offset)/radius);
        //float weight = exp(-0.5*dXY2) * pow(1.0 + dI2/power, - power);
        vec2 krnl = Kernel(vec2(X,Y) - offset);
        float weight = krnl.x * krnl.y * pow(1.0 + dI2/power, - power);

        vec4 lr = LowRes(X,Y);
        lr.yz -= textureLod(CHROMA_raw, ddxddy*(pos+vec2(X,Y)+0.5), 0.0).xy;
        diff += weight * lr;
        weightSum += weight;
    }
    diff /= weightSum;
    c0.yz -= strength * diff.yz;

    return c0;
}
