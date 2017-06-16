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
testSimilarity descript ε a = QC.testProperty (descript++", with ε="++show ε++"")
                                $ epsiloned ε a

