#!/usr/bin/env stack
  {- stack
    script
    --resolver lts-18.16
    --package base
    --package Chart
  -}
import           Graphics.Rendering.Chart.Easy
import           Graphics.Rendering.Chart.Backend.Cairo

signal :: [Double] -> [(Double, Double)]
signal xs =
  [ (x, (sin (x * 3.14159 / 45) + 1) / 2 * (sin (x * 3.14159 / 5))) | x <- xs ]

main = toFile def "example1_big.png" $ do
  layout_title .= "Amplitude Modulation"
  setColors [opaque blue, opaque red]
  plot (line "am" [signal [0, (0.5) .. 400]])
  plot (points "am points" (signal [0, 7 .. 400]))



case dkjf of
  dkj -> value

          otherwise -> other
