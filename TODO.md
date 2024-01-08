# Wish list

- Consider using properties instead of getters/setters.

- Make `vertices` return a vector.

- Simplify building of polygons because the `convex` and `direct` fields are
  not needed by the code for forging mask.

- Additional type for complex mask consisting in the combination of elementary
  apertures and obscurations.

- Affine geometric transforms to geometric objects.

- Broadcast `.*` as `*` for multiplication of geometric objects by a scalar.

- Export `grow` and `shrink`?

- Implement `rotate(obj; by::Angle, around::Point)`.

- Implement `origin(::Type{<:GeometricObject{T}}) where {T} = Point(zero(T), zero(T))`.

- Revise the logic of adding or subtracting a point to a geometric object. If
  the geometric object is seen as a set of points (possibly empty for a box),
  then adding or subtracting a point to a geometric object should yield the set
  corresponding to the point-wise operation between the input point and all
  points of the input set. This may be restricted to `.+` and `.-` operators.

- Change the 50% antialiasing rule to apply for `antialiasing ≤ 0`, not for
  `antialiasing ≤ 1`.

- Check optimality of the forging mask algorithm regarding the use of bounding
  boxes.

- Add constructor for regular polygons. For example, with keyword-based
  constructor: `Polygon(; orient=, radius=, center=, edges=)`.

- Add elliptical shapes.

- Add polyline shapes = connected segments with thickness, cap style, and join
  style.

- Add other cap styles than *butt* for segments (*projecting* and *round*).

- Consider interfacing `Cairo` for drawing masks.

- A segmented pupil may have different piston and/or tip-tilt per segment.
