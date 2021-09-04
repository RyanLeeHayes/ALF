#! /bin/bash

echo "E coli"
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==0) {print $0}}' Corrected.txt

echo "E coli background"
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==1) {print $0}}' Corrected.txt

echo "Most stable E coli mutant"
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==1 && $16<M) {M=$16; MS=$0}}END{print MS}' Corrected.txt

echo "Least stable E coli mutant"
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==1 && $16>M) {M=$16; MS=$0}}END{print MS}' Corrected.txt

echo "AncCCons"
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==15) {print $0}}' Corrected.txt

echo "AncCCons background"
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==14) {print $0}}' Corrected.txt

echo "Most stable AncCCons mutant"
awk 'BEGIN{M=100} {if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==14 && $16<M) {M=$16; MS=$0}}END{print MS}' Corrected.txt

echo "Least stable AncCCons mutant"
awk 'BEGIN{M=-100} {if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15==14 && $16>M) {M=$16; MS=$0}}END{print MS}' Corrected.txt

echo "Most stable chimera"
awk '{if ($16<M) {M=$16; MS=$0}}END{print MS}' Corrected.txt

echo "Least stable chimera"
awk '{if ($16>M) {M=$16; MS=$0}}END{print MS}' Corrected.txt

