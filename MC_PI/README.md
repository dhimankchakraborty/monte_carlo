# Calculating the value of $\pi$ using (Classical) Monte Carlo Method

## 📌 Overview

This project demonstrates how to estimate the value of π using a Monte Carlo method. The idea is based on geometric probability by randomly sampling points inside a unit square.

## 🧠 Mathematical Idea

We consider a unit square and a quarter circle of radius 1 inside it.

The probability that a randomly chosen point lies inside the quarter circle is:

$P = \frac{\text{Area of Quarter Circle}}{\text{Area of Square}} = \frac{\pi/4}{1} = \frac{\pi}{4}$

Thus, we estimate:

$\pi \approx 4 \times \frac{\text{Number of points inside circle}}{\text{Total number of points}}$

## ⚙️ Algorithm

1. Take the number of random points (N) as input
2. Initialize a counter for points inside the circle
3. Repeat (N) times:
    - Generate random numbers ($x, y \in [0,1]$)
    - Check if (x^2 + y^2 $\leq$ 1)
    - If true, increment the counter
4. Estimate $\pi$ using:
$\pi \approx 4 \times \frac{\text{inside count}}{N}$
5. Compute percentage error
