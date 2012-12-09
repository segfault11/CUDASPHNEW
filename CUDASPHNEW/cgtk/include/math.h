/*
 * Utilities for Computer Graphics Applications
 * Copyright (C) 2012 Arno Luebke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef MATH_H
#define MATH_H

typedef float CGTKVec3f[3];

void cgtkVec3fDotProduct(CGTKVec3f u, CGTKVec3f v, float* result );
void cgtkVec3fCrossProduct(CGTKVec3f u, CGTKVec3f v, CGTKVec3f w);
void cgtkVec3fNormalize(CGTKVec3f v);

#endif /* end of include guard: MATH_H */
