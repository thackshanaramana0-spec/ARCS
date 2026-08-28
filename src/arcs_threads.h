#pragma once
// ── One place that decides how many threads ARCS may use ─────────────────────
//
// Before this, std::thread::hardware_concurrency() was called at sixteen
// independent sites (assembly growth, the fallback placement, SA construction,
// the APSP walk, repeat elimination, quality blocks, name coding, the decoder,
// and more). Each of them took every core on the machine, and there was no way
// to ask for less.
//
// That is a problem for anyone running ARCS next to their own work: sixteen
// places each spawning hardware_concurrency() threads will starve whatever else
// is on the box, and on a shared machine it also makes ARCS's own timings
// unreproducible -- two ARCS runs on the same 12 cores contend with each other
// and both report inflated wall times.
//
// ARCS_THREADS=N caps every one of those sites. Unset, the behaviour is exactly
// what it was: hardware_concurrency(). So the default path is unchanged, and
// nothing about the archive depends on this value on any code path that is
// deterministic to begin with.
//
// NOTE, deliberately: the growth loop's output DOES depend on its thread count,
// because threads claim reads by atomic compare-exchange and a contested read
// goes to whoever gets there first (see build_vodbg_pg). Lowering ARCS_THREADS
// therefore changes the pseudogenome slightly, exactly as ARCS_DETERMINISTIC=1
// does -- every result stays lossless, but archives are not comparable across
// different thread counts. Fix the value when benchmarking.
#include <thread>
#include <cstdlib>
#include <algorithm>

inline int arcs_threads() {
    static const int v = [] {
        int hw = (int)std::thread::hardware_concurrency();
        if (hw < 1) hw = 4;                       // hardware_concurrency may report 0
        if (const char* e = std::getenv("ARCS_THREADS")) {
            const int n = std::atoi(e);
            if (n >= 1) return std::min(n, hw);   // never more than the machine has
        }
        return hw;
    }();
    return v;
}
