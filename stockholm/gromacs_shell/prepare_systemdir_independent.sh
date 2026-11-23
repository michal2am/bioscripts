#!/usr/bin/env bash
# Prepare GROMACS run directories from CHARMM-GUI archives named *_sysN.tgz
# Usage: run in the directory containing the archives:  bash setup_charmm_systems.sh
set -euo pipefail
shopt -s nullglob dotglob

archives=( *_sys[0-9]*.tgz )
(( ${#archives[@]} )) || { echo "No *_sysN.tgz archives found in $(pwd)" >&2; exit 1; }

mkdir -p mdp

for tgz in "${archives[@]}"; do
    n=${tgz##*_sys}; n=${n%.tgz}
    [[ $n =~ ^[0-9]+$ ]] || { echo "Skipping $tgz (cannot parse sys number)" >&2; continue; }
    d="sys${n}"
    echo ">> $tgz -> $d/"

    mkdir -p "$d/full_system"
    cp "$tgz" "$d/"
    tar -xzf "$d/$tgz" -C "$d/full_system"

    # CHARMM-GUI packs everything into one top dir (charmm-gui-XXXXXXXXX/) - flatten it
    if [[ ! -d "$d/full_system/gromacs" ]]; then
        top=( "$d"/full_system/*/ )
        if (( ${#top[@]} == 1 )) && [[ -d "${top[0]}gromacs" ]]; then
            mv "${top[0]}"* "$d/full_system/"
            rmdir "${top[0]}"
        fi
    fi

    g="$d/full_system/gromacs"
    [[ -d "$g" ]] || { echo "ERROR: no gromacs/ directory found in $tgz" >&2; exit 1; }

    cp -r "$g/toppar" "$d/"
    files=( "$g"/*.top "$g"/*.ndx "$g"/step5* "$g"/*.itp )
    (( ${#files[@]} )) && cp "${files[@]}" "$d/"

    rm -f "$tgz"
done


echo "Done. Created: $(ls -d sys*/ mdp/ | tr '\n' ' ')"
