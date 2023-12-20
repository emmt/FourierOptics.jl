# The code in this file is patterned after C++ code by Dan Sunday and subjects to
# the following copyright notice.
#
#     Copyright 2000 softSurfer, 2012 Dan Sunday
#     This code may be freely used and modified for any purpose
#     providing that this copyright notice is included with it.
#     SoftSurfer makes no warranty for this code, and cannot be held
#     liable for any real or imagined damage resulting from its use.
#     Users of this code must verify correctness for their application.
#
# Article and original C++ code are available at:
# https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html

"""
    crossing_number_test(P::Point, V::Points) -> bool

yields whether point `P` is inside the polygon defined by vertices `V`
according to the *crossing number* method by Franklin (2000).

"""
crossing_number_test(P::Point, V::Points) = isodd(crossing_number(P, V))

function crossing_number(P::Point, V::Points)
    cn = 0 # initialize crossing number counter
    A = last(V) # initialize A, the 1st point of edges in cyclic list of points
    @inbounds for B in V # loop over B, the 2nd point of edges
        if  A.y ≤ P.y < B.y || # an upward crossing
            B.y ≤ P.y < A.y    # a downward crossing

            # Compute the actual edge-ray intersect x-coordinate.
            t = (P.y  - A.y)/(B.y - A.y)
            if P.x < A.x + t*(B.x - A.x) # P.x < intersect
                cn += 1 # a valid crossing of y=P.y right of P.x
            end
        end
        A = B # update 1st point of next edge
    end
    return cn # outside if even, inside if odd
end

"""
    winding_number_test(P::Point, V::Points) -> bool

yields whether point `P` is inside the polygon defined by vertices `V`
according to the *winding number* method by Dan Sunday ("Inclusion of a Point
in a Polygon", 2001).

"""
winding_number_test(P::Point, V::Points) = !iszero(winding_number(P, V))

function winding_number(P::Point, V::Points)
    wn = 0 # initialize winding number counter
    A = last(V) # initialize A, the 1st point of edges in cyclic list of points
    @inbounds for B in V # loop over B, the 2nd point of edges
        if A.y ≤ P.y
            if B.y > P.y  # an upward crossing
                s = cross(A, B, P)
                if s > zero(s) # P left of edge
                    wn += 1 # have a valid up intersect
                end
            end
        else # A.y > P.y (no test needed unless NaN)
            if B.y ≤ P.y # a downward crossing
                s = cross(A, B, P)
                if s < zero(s) # P right of  edge
                    wn -= 1 # have a valid down intersect
                end
            end
        end
        A = B # update 1st point of next edge
    end
    return wn # 0 only when P is outside
end

"""
    cross(A::Point, B::Point, C::Point) -> val

yields the cross vectorial product of `AB` by `AC`.

The returned value can be used to determine the position of point `C`
relatively to the infinite line defined by `(A,B)`:

- if `val > 0`, then `C` is left of the line through `A` and `B`;
- if `val = 0`, then `C` is on the line through `A` and `B`;
- if `val < 0`, then `C` is right of the line through `A` and `B`.

See: Algorithm 1 "Area of Triangles and Polygons"

"""
cross(A::Point, B::Point, C::Point) =
    (B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y)
