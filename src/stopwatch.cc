/*==============================================================================
  BSD 2-Clause License

  Copyright (c) 2025 Shogo OKADA (shogo.okada@kek.jp)
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
==============================================================================*/
/*============================================================================
Copyright 2017-2022 Koichi Murakami

Distributed under the OSI-approved BSD License (the "License");
see accompanying file LICENSE for details.

This software is distributed WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the License for more information.
============================================================================*/
#include "stopwatch.hh"

// --------------------------------------------------------------------------
Stopwatch::Stopwatch()
{
  Reset();
}

// --------------------------------------------------------------------------
void Stopwatch::Reset()
{
  start_clock_ = times(&start_time_);
}

// --------------------------------------------------------------------------
void Stopwatch::Split()
{
  end_clock_ = times(&end_time_);
}

// --------------------------------------------------------------------------
double Stopwatch::GetRealElapsed() const
{
  double diff = end_clock_ - start_clock_;
  return diff / sysconf(_SC_CLK_TCK);
}

// --------------------------------------------------------------------------
double Stopwatch::GetSystemElapsed() const
{
  double diff = end_time_.tms_stime - start_time_.tms_stime;
  return diff / sysconf(_SC_CLK_TCK);
}

// --------------------------------------------------------------------------
double Stopwatch::GetUserElapsed() const
{
  double diff = end_time_.tms_utime - start_time_.tms_utime;
  return diff / sysconf(_SC_CLK_TCK);
}

// --------------------------------------------------------------------------
const char* Stopwatch::GetClockTime() const
{
   time_t timer;
   time(&timer);

   return ctime(&timer);
}
