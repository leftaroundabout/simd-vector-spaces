-- |
-- Module      : Main
-- Copyright   : (c) Justus Sagemüller 2017
-- License     : GPL v3
-- 
-- Maintainer  : (@) sagemueller $ geo.uni-koeln.de
-- Stability   : experimental
-- Portability : portable
-- 
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE LambdaCase       #-}

module Main where

import Math.DivisionAlgebra.SIMD
import Data.VectorSpace.SIMD

import Data.AdditiveGroup
import Data.VectorSpace

import Test.Tasty
import Test.Tasty.QuickCheck as QC


main = defaultMain tests



prop_distrib :: CheckableVect v => Scalar v -> v -> v -> Similarity
prop_distrib μ v w = μ *^ v ^+^ μ *^ w ≈ μ *^ (v^+^w)


tests :: TestTree
tests = testGroup "Vector-space identities"
  [testGroup "(checked by QuickCheck)"
     [ testSimilarity "distrib @(ℝ,ℝ)" 1e-14
        (prop_distrib :: ℝ -> (ℝ,ℝ) -> (ℝ,ℝ) -> Similarity)
     , testSimilarity "distrib @ℝ³" 1e-14
        (prop_distrib :: ℝ -> ℝ³ -> ℝ³ -> Similarity)
     ]
  ]



type CheckableVect v = (InnerSpace v, Num (Scalar v), Ord (Scalar v))
newtype Similarity = Similarity {evalSimilarity :: ℝ -> Bool}

infix 4 ≈
(≈) :: CheckableVect v => v -> v -> Similarity
v ≈ w = Similarity $ \ε -> magnitudeSq (v^-^w) * fromInteger (round $ recip ε^2)
                         <= (magnitudeSq v + magnitudeSq w)


class (QC.Testable (AfterEpsilon t)) => SimTestable t where
  type AfterEpsilon t :: *
  epsiloned :: ℝ -> t -> AfterEpsilon t
instance SimTestable Similarity where
  type AfterEpsilon Similarity = Bool
  epsiloned ε (Similarity s) = s ε
instance (QC.Arbitrary a, Show a, SimTestable sim) => SimTestable (a -> sim) where
  type AfterEpsilon (a -> sim) = a -> AfterEpsilon sim
  epsiloned ε sim = epsiloned ε . sim

testSimilarity :: SimTestable a => TestName -> ℝ -> a -> TestTree
testSimilarity descript ε a = testGroup descript
       [ QC.testProperty ("with ε="++showOOM ε++"") $ epsiloned ε a
       , QC.testProperty ("with ε="++showOOM (ε/100)++"") . QC.expectFailure
               $ epsiloned (ε/100) a ]



showOOM :: Double -> String
showOOM n = case m₁₀ of
   1 -> "10"++showSup e₁₀
   _ -> show m₁₀ ++ "×10" ++ showSup e₁₀
 where e₁₀ = floor $ logBase 10 n
       m₁₀ = round $ n / 10^^e₁₀
       showSup = map (\case
        {'0'->'⁰';'1'->'¹';'2'->'²';'3'->'³';'4'->'⁴';'5'->'⁵';'6'->'⁶'
        ;'7'->'⁷';'8'->'⁸';'9'->'⁹';'-'->'⁻';c->c})
        . show
