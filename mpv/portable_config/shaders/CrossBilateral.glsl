// CrossBilateral by Shiandow
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


//!HOOK CHROMA
//!BIND HOOKED
//!BIND LUMA
//!SAVE DOWNSCALEDLUMAX
//!HEIGHT LUMA.h
//!WHEN CHROMA.w LUMA.w <
//!COMPONENTS 4

// -- Downscaling --
#define offset   (-vec2(0.0, 0.0)*LUMA_size*CHROMA_pt) // -chromaOffset*LUMA_size/CHROMA_size

#define dxdy     (vec2(CHROMA_pt.x, LUMA_pt.y)) // 1/output_size
#define ddxddy   (LUMA_pt)                      // 1/input_size

#define factor   ((ddxddy*vec2(CHROMA_size.x, LUMA_size.y))[axis])

#define axis 0

#define Kernel(x) clamp(0.5 + (0.5 - abs(x)) / factor, 0.0, 1.0)
#define taps (1.0 + factor)

vec4 hook() {
    // Calculate bounds
    float low  = floor((LUMA_pos - 0.5*taps*dxdy) * LUMA_size - offset + 0.5)[axis];
    float high = floor((LUMA_pos + 0.5*taps*dxdy) * LUMA_size - offset + 0.5)[axis];

    float W = 0.0;
    vec2 avg = vec2(0);
    vec2 pos = LUMA_pos;

    for (float k = 0.0; k < high - low; k++) {
        pos[axis] = ddxddy[axis] * (k + low + 0.5);
        float rel = (pos[axis] - LUMA_pos[axis])*vec2(CHROMA_size.x, LUMA_size.y)[axis] + offset[axis]*factor;
        float w = Kernel(rel);

        avg += w * vec2(textureLod(LUMA_raw, pos, 0.0).x, pow(textureLod(LUMA_raw, pos, 0.0).x, 2.0));
        W += w;
    }
    avg /= vec2(W);

    return vec4(avg[0], avg[1] - avg[0]*avg[0], 0, 0);
}

//!HOOK CHROMA
//!BIND HOOKED
//!BIND DOWNSCALEDLUMAX
//!SAVE LOWRES_YUV
//!COMPONENTS 4

// -- Downscaling --
#define offset   (-vec2(0.0, 0.0)*DOWNSCALEDLUMAX_size*CHROMA_pt)

#define dxdy     (CHROMA_pt)
#define ddxddy   (DOWNSCALEDLUMAX_pt)

#define factor   ((ddxddy*CHROMA_size)[axis])

#define axis 1

#define Kernel(x) clamp(0.5 + (0.5 - abs(x)) / factor, 0.0, 1.0)
#define taps (1.0 + factor)

vec4 hook() {
    // Calculate bounds
    float low  = floor((DOWNSCALEDLUMAX_pos - 0.5*taps*dxdy) * DOWNSCALEDLUMAX_size - offset + 0.5)[axis];
    float high = floor((DOWNSCALEDLUMAX_pos + 0.5*taps*dxdy) * DOWNSCALEDLUMAX_size - offset + 0.5)[axis];

    float W = 0.0;
    vec2 avg = vec2(0);
    vec2 pos = DOWNSCALEDLUMAX_pos;

    for (float k = 0.0; k < high - low; k++) {
        pos[axis] = ddxddy[axis] * (k + low + 0.5);
        float rel = (pos[axis] - DOWNSCALEDLUMAX_pos[axis])*CHROMA_size[axis] + offset[axis]*factor;
        float w = Kernel(rel);

        avg += w * vec2(textureLod(DOWNSCALEDLUMAX_raw, pos, 0.0).x, textureLod(DOWNSCALEDLUMAX_raw, pos, 0.0).y + pow(textureLod(DOWNSCALEDLUMAX_raw, pos, 0.0).x, 2.0));
        W += w;
    }
    avg /= vec2(W);

    return vec4(avg[0], texture(CHROMA_raw, CHROMA_pos).xy, avg[1]-avg[0]*avg[0]);
}

//!HOOK CHROMA
//!BIND HOOKED
//!BIND LUMA
//!BIND LOWRES_YUV
//!WIDTH LUMA.w
//!HEIGHT LUMA.h

// -- CrossBilateral --

// -- Misc --
#define power  0.5

// -- Convenience --
#define sqr(x) dot(x,x)
#define noise  0.1
#define bitnoise 1.0/(2.0*255.0)
#define radius 1.0
#define chromaOffset vec2(0.0, 0.0)

// -- Window Size --
#define taps 3
#define even (float(taps) - 2.0 * floor(float(taps) / 2.0) == 0.0)
#define minX int(1.0-ceil(float(taps)/2.0))
#define maxX int(floor(float(taps)/2.0))

// -- Input processing --
// Luma value
//#define GetLuma(x,y)   LOWRES_YUV_tex(HOOKED_pos + NATIVE_pt*vec2(x,y))[0]
// Chroma value
#define GetChroma(x,y) LOWRES_YUV_tex(LOWRES_YUV_pt*(pos+vec2(x,y)+vec2(0.5)))

#define localVar 0.0//sqr(bitnoise)
#define fixLuma 0.0

#define C(i,j) (inversesqrt(sqr(noise) + X[i].w + X[j].w) * exp(-0.5*(sqr(X[i].x - X[j].x)/(sqr(noise) + X[i].w + X[j].w) + sqr((coords[i] - coords[j])/radius))) + fixLuma * (X[i].x - y) * (X[j].x - y))
#define c(i) (inversesqrt(sqr(noise) + X[i].w + localVar) * exp(-0.5*(sqr(X[i].x - y)/(sqr(noise) + X[i].w + localVar) + sqr((coords[i] - offset)/radius))))

#define N (taps*taps - 1)

#define f1(i) {b[i] = c(i) - c(N) - C(i,N) + C(N,N); \
        for (int j=i; j<N; j++) M[(i)*N + j] = C(i,j) - C(i,N) - C(j,N) + C(N,N);}

#define f2(i) for (int j=i+1; j<N; j++) {b[j] -= b[i] * M[(i)*N + j] / M[(i)*N + (i)]; \
        for (int k=j; k<N; k++) M[j*N + k] -= M[(i)*N + k] * M[(i)*N + j] / M[(i)*N + (i)];}

#define I2(f, n) f(n) f(n+1)
#define I4(f, n) I2(f, n) I2(f, n+2)
#define I8(f, n) I4(f, n) I4(f, n+4)

vec4 hook() {
    float y = LUMA_tex(LUMA_pos).x;

    // Calculate position
    vec2 pos = LUMA_pos * LOWRES_YUV_size - chromaOffset - vec2(0.5);
    vec2 offset = pos - (even ? floor(pos) : round(pos));
    pos -= offset;

    vec2 coords[N+1];
    vec4 X[N+1];
    int i=0;
    for (int xx = minX; xx <= maxX; xx++)
    for (int yy = minX; yy <= maxX; yy++) 
        if (!(xx == 0 && yy == 0)) {
            coords[i] = vec2(xx,yy);
            X[i++] = GetChroma(xx, yy);
        }

    coords[N] = vec2(0,0);
    X[N] = GetChroma(0,0);

    float M[N*N];
    float b[N];

    I8(f1, 0)

    I8(f2, 0)

    //float w[N];
    //float det = 1;
    //float Tr = 0;
    vec4 interp = X[N];
    for (i=N-1; i>=0; i--) {
        //w[i] = b[i];
        for (int j=i+1; j<N; j++) {
            b[i] -= M[i*N + j] * b[j];
        }
        b[i] /= M[i*N + i];
        interp += b[i] * (X[i] - X[N]);
        //det *= M[i*N + i];
        //Tr += M[i*N + i];
    }

    return interp.yzxx;
}
