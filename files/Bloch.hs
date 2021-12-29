#!/bin/env stack
{- stack
  script
  --resolver lts-18.19
  --package hmatrix
  --package lens
  --package mtl
  --package aeson
  --package hvega
-}

{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}

module Mag where

import           Control.Lens            hiding ( transform )
import           Control.Monad.State.Lazy
import           Data.Aeson                     ( encodeFile )
import           Data.Bifunctor
import           Data.List.NonEmpty             ( NonEmpty(..) )
import qualified Data.List.NonEmpty            as Lne
import           Graphics.Vega.VegaLite
import           Numeric.LinearAlgebra
import           Prelude                 hiding ( (<>) )

type Mag = Vector R
type T1_s = Double
type T2_s = Double
type Time_s = Double
type Freq_Hz = Double
type Angle_rad = Double
type Axis = [Double]

type BlochS = (Time_s, Mag)

type BlochST = State (NonEmpty BlochS) ()

initBlochS :: Time_s -> [Double] -> BlochS
initBlochS t es = (t, fromList es)

updateBlochST :: (BlochS -> BlochS) -> BlochST
updateBlochST f = do
  (s :| _) <- get
  modify (Lne.cons (f s))

runBlochST :: BlochST -> BlochS -> ([Time_s], [Mag])
runBlochST st s0 =
  unzip $ (s0 :) $ Lne.toList $ Lne.reverse $ execState st (s0 :| [])

blockST
  :: Axis -> Angle_rad -> Freq_Hz -> T1_s -> T2_s -> Time_s -> Int -> BlochST
blockST ax th freq t1 t2 dt n = do
  instPulseST ax th
  fidST dt freq t1 t2 n

instPulseST :: Axis -> Angle_rad -> BlochST
instPulseST ax th = updateBlochST (second $ instPulseM (fromList ax) th)

instPulseM :: Vector R -> Angle_rad -> Mag -> Mag
instPulseM u th m = rot u th #> m

fidST :: Time_s -> Freq_Hz -> T1_s -> T2_s -> Int -> BlochST
fidST dt freq t1 t2 n = do
  let (precs, decay, b) = freeprecess dt freq t1 t2
      a                 = precs <> decay
  replicateM_ n $ updateBlochST (bimap (+ dt) (\m -> a #> m + b))

fidM :: Time_s -> Freq_Hz -> T1_s -> T2_s -> Mag -> Mag
fidM dt freq t1 t2 m =
  let (precs, decay, b) = freeprecess dt freq t1 t2
      a                 = precs <> decay
  in  a #> m + b

freeprecessM :: Time_s -> Freq_Hz -> T1_s -> T2_s -> Mag -> Mag
freeprecessM dt freq t1 t2 m = (precs <> decay) #> m + recov
  where (precs, decay, recov) = freeprecess dt freq t1 t2

freeprecess
  :: Time_s -> Freq_Hz -> T1_s -> T2_s -> (Matrix R, Matrix R, Vector R)
freeprecess dt freq t1 t2 = (precs, decay, recov)
 where
  e1    = exp $ -dt / t1
  e2    = exp $ -dt / t2
  decay = diagl [e2, e2, e1]
  precs = precess dt [0, 0, 1] freq
  recov = fromList [0, 0, 1 - e1]

precess :: Time_s -> Axis -> Freq_Hz -> Matrix R
precess dt ax freq = rot (fromList ax) (2 * pi * freq * dt)

rot :: Vector R -> Angle_rad -> Matrix R
rot u th = scale c (ident 3) + scale (sin th) (crossProdMat u) + scale
  (1 - c)
  (outer u u)
  where c = cos th

crossProdMat :: Vector R -> Matrix R
crossProdMat v = fromColumns $ map (cross v) (toColumns $ ident 3)

plotComponents :: FilePath -> [PropertySpec] -> ([Time_s], [Mag]) -> IO ()
plotComponents fname props xs = do
  let vega = toVegaLite $ vlDataToPropSpecs (toVLData xs) ++ props
  encodeFile fname $ fromVL vega

toVLData :: ([Time_s], [Mag]) -> Data
toVLData (ts, ms) =
  dataFromColumns []
    . dataColumn "time" (Numbers ts)
    . dataColumn "M_x"  (Numbers mx)
    . dataColumn "M_y"  (Numbers my)
    . dataColumn "M_z"  (Numbers mz)
    $ []
 where
  mx = map (! 0) ms
  my = map (! 1) ms
  mz = map (! 2) ms

vlDataToPropSpecs :: Data -> [PropertySpec]
vlDataToPropSpecs dat =
  [ dat
  , trans []
  , mark Line [MTooltip TTEncoding, MClip True]
  , enc []
  , width 800
  , height 450
  ]
 where
  trans = transform . fold ["M_x", "M_y", "M_z"]
  enc =
    encoding
      . position
          X
          [ PName "time"
          , PTitle "Time [s]"
          , PmType Quantitative
          -- , PScale [SDomain (DNumbers [0, 1.0])]
          ]
      . position
          Y
          [PName "value", PTitle "Magnetization [a.u.]", PmType Quantitative]
      . color [MName "key", MTitle "Component"]

main_fid :: IO ()
main_fid = do
  let dt     = 1e-3
      m0     = [1, 0, 0]
      freq   = 10
      t1     = 600e-3
      t2     = 100e-3
      n      = 1000
      res = runBlochST (fidST dt freq t1 t2 n) (initBlochS 0 m0)
  plotComponents "fid.vg.json" [title "Free induction decay" []] res

main_satrecov :: IO ()
main_satrecov = do
  let nTR       = 10
      ey        = [0, 1, 0]
      flipAngle = pi / 3
      freq      = 0
      t1        = 600e-3
      t2        = 100e-3
      dt        = 1e-3
      nFid      = 500
      m0        = [0, 0, 1]
      res       = runBlochST
        (replicateM_ nTR $ blockST ey flipAngle freq t1 t2 dt nFid)
        (initBlochS 0 m0)
  plotComponents "satrecov.vg.json" [title "Saturation recovery" []] res

main_se :: IO ()
main_se = do
  let ey   = [0, 1, 0]
      ex   = [1, 0, 0]
      freq = 10
      t1   = 600e-3
      t2   = 100e-3
      dt   = 1e-3
      te   = 50e-3
      tR   = 500e-3
      nTR  = 1
      m0   = [0, 0, 1]
      res  = runBlochST
        (do
          replicateM_ nTR $ do
            blockST ey (pi / 2) freq t1 t2 dt (floor $ te / 2 / dt)
            blockST ex pi       freq t1 t2 dt (floor $ (tR - te / 2) / dt)
        )
        (initBlochS 0 m0)
  plotComponents "se.vg.json" [title "Spin echo" []] res

main :: IO ()
main = do
  main_fid
  main_satrecov
  main_se
