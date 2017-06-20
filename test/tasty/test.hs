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
import Data.VectorSpace.SIMD.FiniteSupportedSequence

import qualified Data.Primitive.SIMD as SIMD

import Data.AdditiveGroup
import Data.VectorSpace

import Test.Tasty
import Test.Tasty.QuickCheck as QC

import GHC.Exts (IsList (..))


main = defaultMain tests



prop_distrib :: CheckableVect v => Scalar v -> v -> v -> Similarity
prop_distrib μ v w = μ *^ v ^+^ μ *^ w ≈ μ *^ (v^+^w)

prop_listScalarProd :: ( InnerSpace v, Num (Scalar v), CheckableVect (Scalar v)
                   , IsList v, Item v ~ Scalar v )
                       => v -> v -> Similarity
prop_listScalarProd v w = v<.>w ≈ sum (zipWith (*) (toList v) (toList w))


tests :: TestTree
tests = testGroup "Vector-space identities"
  [testGroup "(checked by QuickCheck)"
     [ testSimilarity "distrib @(ℝ,ℝ)" 1e-14
        (prop_distrib :: ℝ -> (ℝ,ℝ) -> (ℝ,ℝ) -> Similarity)
     , testSimilarity "distrib @ℝ³" 1e-14
        (prop_distrib :: ℝ -> ℝ³ -> ℝ³ -> Similarity)
     , testSimilarity "distrib @ℝⁿ" 1e-14
        (prop_distrib :: ℝ -> FinSuppSeq ℝ ℝ -> FinSuppSeq ℝ ℝ -> Similarity)
     , testSimilarity "scalarProd @ℝⁿ" 0
        (prop_listScalarProd :: FinSuppSeq ℝ ℝ -> FinSuppSeq ℝ ℝ -> Similarity)
     , testSimilarity "distrib @ℝ⁴ⁿ" 1e-14
        (prop_distrib :: ℝ -> FinSuppSeq ℝX4 ℝ -> FinSuppSeq ℝX4 ℝ -> Similarity)
     , testSimilarity "scalarProd @ℝ⁴ⁿ" 1e-14
        (prop_listScalarProd :: FinSuppSeq ℝX4 ℝ -> FinSuppSeq ℝX4 ℝ -> Similarity)
     ]
  ]



type CheckableVect v = (InnerSpace v, Num (Scalar v), Ord (Scalar v))
newtype Similarity = Similarity {evalSimilarity :: ℝ -> Bool}

infix 4 ≈
(≈) :: CheckableVect v => v -> v -> Similarity
v ≈ w = Similarity simvw
 where simvw 0 = magnitudeSq (v^-^w) == 0
       simvw ε = magnitudeSq (v^-^w) * fromInteger (round $ recip ε^2)
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
      $ QC.testProperty ("with ε="++showOOM ε++"") (epsiloned ε a)
      : [ QC.testProperty ("with ε="++showOOM (ε/100)++"") . QC.expectFailure
               $ epsiloned (ε/100) a
        | ε>0 ]



showOOM :: Double -> String
showOOM 0 = "0"
showOOM n = case m₁₀ of
   1 -> "10"++showSup e₁₀
   _ -> show m₁₀ ++ "×10" ++ showSup e₁₀
 where e₁₀ = floor $ logBase 10 n
       m₁₀ = round $ n / 10^^e₁₀
       showSup = map (\case
        {'0'->'⁰';'1'->'¹';'2'->'²';'3'->'³';'4'->'⁴';'5'->'⁵';'6'->'⁶'
        ;'7'->'⁷';'8'->'⁸';'9'->'⁹';'-'->'⁻';c->c})
        . show


type ℝX4 = SIMD.DoubleX4
