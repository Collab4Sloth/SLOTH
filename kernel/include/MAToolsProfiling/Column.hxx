/**
 * @file Column.hxx
 * @author Raphaël Prat (raphael.prat@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-09-08
 *
 * Copyright CEA (C) 2025
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#pragma once
#include <iostream>
#include <string>

// variables
const int cWidth = 20;   // Column width
const int nColumns = 6;  // Number of columns
const std::string cName[nColumns] = {"number Of calls", "min (s)", "mean (s)", "max (s)",
                                     "time ratio (%)",  "imb (%)"};  // [1-Imax/Imean]%
