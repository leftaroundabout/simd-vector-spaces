-- |
-- Module      : Data.VectorSpace.SIMD
-- Copyright   : (c) Justus Sagemüller 2017
-- License     : GPL v3
-- 
-- Maintainer  : (@) sagemueller $ geo.uni-koeln.de
-- Stability   : experimental
-- Portability : portable
-- 
{-# LANGUAGE TypeFamilies            #-}
{-# LANGUAGE TypeOperators           #-}
{-# LANGUAGE FlexibleInstances       #-}
{-# LANGUAGE FlexibleContexts        #-}
{-# LANGUAGE MultiParamTypeClasses   #-}
{-# LANGUAGE DefaultSignatures       #-}
{-# LANGUAGE ScopedTypeVariables     #-}
{-# LANGUAGE UnicodeSyntax           #-}
{-# LANGUAGE CPP                     #-}

module Data.VectorSpace.SIMD (Complex) where

import Math.DivisionAlgebra.SIMD

import Data.AdditiveGroup
import Data.VectorSpace

import qualified Data.Primitive.SIMD as SIMD


type ℝ = Double
newtype ℝ² = ℝ² SIMD.DoubleX2

#define SIMDVSInstances(V,F)       \
instance AdditiveGroup (V) where {  \
  V a ^+^ V b = V (a+b);             \
  zeroV = V 0;                        \
  negateV (V a) = V $ negate a };      \
instance VectorSpace V where {          \
  type Scalar (V) = (F);                 \
  μ *^ V a = V $ SIMD.broadcastVector μ * a }

SIMDVSInstances(ℝ², ℝ)


data DirectSum v w = DirectSum !v !w

infixl 6 ⊕
type family (⊕) v w where
  ℝ ⊕ ℝ = ℝ²
  v ⊕ w = DirectSum v w


class KnownDirectSum v w where
  (⊕) :: v -> w -> v⊕w
  default (⊕) :: v -> w -> DirectSum v w
  (⊕) = DirectSum
  split :: v⊕w -> (v,w)
  default split :: DirectSum v w -> (v,w)
  split (DirectSum v w) = (v,w)

instance KnownDirectSum ℝ ℝ where
  {-# INLINE (⊕) #-}
  (⊕) x y = ℝ² $ SIMD.packVector (x,y)
  {-# INLINE split #-}
  split (ℝ² a) = SIMD.unpackVector a


