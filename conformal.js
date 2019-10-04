// TODO: Apply better math.js idioms file-wide.

function densify(vector, segments){

    numel = math.squeeze(math.size(vector))
    if (!math.isInteger(numel)){
        Error("Cannot densify non-vectors");
    }
    
    const shifted = vector.subset(math.index(math.concat(math.range(1, math.squeeze(math.size(vector))), [0])));
    
    return math.flatten(
        math.add(
            math.multiply(
                math.transpose([vector]),
                math.ones(1, segments)
            ),
            math.multiply( // Hacky implementation of outer product
                math.transpose([math.subtract(shifted, vector)]),
                [math.range(0, segments)],
                1 / segments,
            )
        )
    )
}



function geod_fwd_phi_1(z, z0, z1){
    return math.map(z, function (x) {
        // z -> sqrt((z - z1) / (z - z0))
        w = math.multiply(
            math.evaluate("1i"),
            math.sqrt(
                math.divide(
                    math.subtract(x, z1),
                    math.subtract(x, z0))
            )
        );

        return w;
    })
}

function geod_fwd_phi_i(z, zi){
    b = math.pow(math.abs(zi), 2) / zi.re;
    csq = math.pow(math.pow(math.abs(zi), 2) / zi.im, 2);
    
    return math.map(z, function (x) {
        var w = x;
        // Ignore first automorphism if b is infinite
        if (math.abs(zi.re) > 0){
            // x is infinite, maps to -b
            if (isNaN(x.re) || isNaN(x.im)){ w = math.complex(-b); }
            // z -> z / (1 - z/b)
            else {
                w = math.divide(
                    x,
                    math.subtract(1, math.divide(x, b))
                )
            }
        }
        
        // z -> sqrt(z^2 + c^2)
        w = math.multiply(
            w,
            math.sqrt(
                math.add(
                    1,
                    math.divide(
                        csq,
                        math.pow(w, 2)
                    )
                )
            )
        )

        return w;
    })
    
}

function geod_fwd_phi_np1(z, zn){
    return math.map(z, function (x) {
    
        // z -> z / (1 - z/zn)
        if (isNaN(x.re) || isNaN(x.im)){ w = math.mul(-1, zn); }
        else {
            w = math.divide(
                x,
                math.subtract(1, math.divide(x, zn))
            )
        }
        
        // z -> -z^2
        w = math.multiply(
            -1,
            math.pow(w, 2)
        )
        
        return w;
    })
}


function geod_inv_phi_1(z, z0, z1){
    return math.map(z, function (x) {
        // z -> (z1 + z0*z^2) / (z^2 + 1)
        sq = math.pow(x, 2)
        w = math.divide(
            math.add(
                z1,
                math.multiply(z0, sq)
            ),
            math.add(
                sq,
                1
            )
        )
        
        return w;
    })
}

function geod_inv_phi_i(z, zi){
    b = math.pow(math.abs(zi), 2) / zi.re;
    csq = math.pow(math.pow(math.abs(zi), 2) / zi.im, 2);
    
    return math.map(z, function (x) {

        // z -> z sqrt(1 - c^2 z^-2)
        w = math.multiply(
            x,
            math.sqrt(
                math.subtract(
                    1,
                    math.divide(
                        csq,
                        math.pow(x, 2)
                    )
                )
            )
        )

        // Ignore first automorphism if b is infinite
        if (math.abs(zi.re) > 0){
            if (isNaN(w.re) || isNaN(w.im)){ w = math.complex(b); }
            // z -> z / (1 + z/b)
            else {
                w = math.divide(
                    w,
                    math.add(1, math.divide(w, b))
                )
            }
        }

        return w;
    })
}

function geod_inv_phi_np1(z, zn){

    return math.map(z, function (x) {
        // inverse of z -> -z^2
        w = math.multiply(
            math.evaluate("1i"),
            math.sqrt(x)
        )
    
        // inverse of z -> z / (1 - z/zn)
        if (isNaN(w.re) || isNaN(w.im)){ w = zn; }
        else {
            w = math.divide(
                w,
                math.add(1, math.divide(w, zn))
            )
        }
        
        return w;
    })
}

function geod_fwd(z, zeta){
    n = math.squeeze(zeta.size())
    z = math.evaluate("fn(v, p[1], p[2])", {v: z, p: zeta, fn: geod_fwd_phi_1})
    for (var i=3; i < n; ++i){
        z = math.evaluate("fn(v, p[" + String(i) + "])", {v: z, p: zeta, fn: geod_fwd_phi_i})
    }
    z = math.evaluate("fn(v, p[end])", {v: z, p: zeta, fn: geod_fwd_phi_np1})
   
    return z;
}

function geod_inv(z, zeta){
    n = math.squeeze(zeta.size())
    z = math.evaluate("fn(v, p[end])", {v: z, p: zeta, fn: geod_inv_phi_np1})
    for (var i=n-1; i > 2; --i){
        z = math.evaluate("fn(v, p[" + String(i) + "])", {v: z, p: zeta, fn: geod_inv_phi_i})
    }
    z = math.evaluate("fn(v, p[1], p[2])", {v: z, p: zeta, fn: geod_inv_phi_1})
   
    return z;
}


function udisk2uhp(z, fromZero){
    return math.evaluate("(z * conj(fromZero) - fromZero) ./ (z - 1)", {z: grid, fromZero: fromZero});
}

function uhp2udisk(z, toZero){
    return math.evaluate("(z - toZero) ./ (z - conj(toZero))", {z: grid, toZero: toZero});
}
