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
#ifndef CLOCK_H
#define CLOCK_H

#ifdef __cplusplus
extern "C"
{
#endif

/* 
**	Starts/restarts the clock 
*/
void cgtkClockStart();
/* 
**	Gets the time past since the clock was started/restarted in ms.
*/
double cgtkClockElapsed();
/*
**	Dumps the time past since the clock was started/restarted in ms.
*/
void cgtkClockDumpElapsed();

#ifdef __cplusplus
}
#endif

#endif /*  end of include guard: clock.h */