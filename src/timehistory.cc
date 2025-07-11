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
#include <iostream>
#include <iomanip>
#include <mutex>
#include "timehistory.hh"

// --------------------------------------------------------------------------
namespace {

std::mutex mtx;

} // end of namespace

// --------------------------------------------------------------------------
TimeHistory* TimeHistory::GetTimeHistory()
{
  static TimeHistory thistory;
  return &thistory;
}

// --------------------------------------------------------------------------
TimeHistory::TimeHistory()
  : sw_(), t0_(0.)
{
  split_history_.clear();

  sw_.Split();
  t0_ = sw_.GetRealElapsed();
}

// --------------------------------------------------------------------------
void TimeHistory::TakeSplit(const std::string& key)
{
  sw_.Split();
  double split = sw_.GetRealElapsed();
  ::mtx.lock();
  split_history_[key] = split - t0_;
  ::mtx.unlock();
}

// --------------------------------------------------------------------------
double TimeHistory::TakeSplit()
{
  sw_.Split();
  double split = sw_.GetRealElapsed();
  double t1 = split - t0_;
  return t1;
}

// --------------------------------------------------------------------------
bool TimeHistory::FindAKey(const std::string& key) const
{
  std::map<std::string, double>::const_iterator itr;
  itr = split_history_.find(key);
  if ( itr != split_history_.end() ) return true;
  else return false;
}

// --------------------------------------------------------------------------
double TimeHistory::GetTime(const std::string& key) const
{
  std::map<std::string, double>::const_iterator itr;
  itr = split_history_.find(key);
  if ( itr != split_history_.end() ) {
    return itr-> second;
  } else {
    std::cout << "[WARNING] TimeHistory::GetTime() cannot find a key. "
              << key << std::endl;
    return 0.;
  }
}

// --------------------------------------------------------------------------
void TimeHistory::ShowHistory(const std::string& key) const
{
  ::mtx.lock();
  std::map<std::string, double>::const_iterator itr;
  itr = split_history_.find(key);
  if ( itr != split_history_.end() ) {
    std::cout << "[" << itr-> first << "] : "
              << itr-> second << "s" << std::endl;
  } else {
    std::cout << "[WARNING] TimeHistory::ShowHistory() cannot find a key. "
              << key << std::endl;
  }
  ::mtx.unlock();
}

// --------------------------------------------------------------------------
void TimeHistory::ShowAllHistories() const
{
  ::mtx.lock();
  std::multimap<double, std::string> histories_by_time;
  std::map<std::string, double>::const_iterator itr;
  for ( itr = split_history_.begin(); itr != split_history_.end(); ++itr) {
    histories_by_time.insert(std::make_pair(itr->second, itr->first));
  }

  std::cout << "* All time histories" << std::endl;

  std::multimap<double, std::string>::const_iterator itr2;
  for ( itr2 = histories_by_time.begin();
        itr2 != histories_by_time.end();  ++itr2) {
    std::cout << "[" << itr2-> second << "] : "
              << std::fixed << std::setprecision(3)
              << itr2-> first << " s" << std::endl;
  }
  ::mtx.unlock();
}

// --------------------------------------------------------------------------
void TimeHistory::ShowClock(const std::string& prefix) const
{
  ::mtx.lock();
  std::cout << prefix << " " << sw_.GetClockTime() << std::flush;
  ::mtx.unlock();
}
