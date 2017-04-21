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
newtype ℝ³ = ℝ³ SIMD.DoubleX4
newtype ℝ⁴ = ℝ⁴ SIMD.DoubleX4

#define SIMDVSInstances(V,F)       \
instance AdditiveGroup (V) where {  \
  V a ^+^ V b = V (a+b);             \
  zeroV = V 0;                        \
  negateV (V a) = V $ negate a };      \
instance VectorSpace V where {          \
  type Scalar (V) = (F);                 \
  μ *^ V a = V $ SIMD.broadcastVector μ * a }

SIMDVSInstances(ℝ², ℝ)
SIMDVSInstances(ℝ³, ℝ)
SIMDVSInstances(ℝ⁴, ℝ)


data DirectSum v w = DirectSum !v !w

infixl 6 ⊕
type family (⊕) v w where
  ℝ ⊕ ℝ = ℝ²
  ℝ ⊕ ℝ² = ℝ³
  ℝ² ⊕ ℝ = ℝ³
  ℝ ⊕ ℝ³ = ℝ⁴
  ℝ³ ⊕ ℝ = ℝ⁴
  ℝ² ⊕ ℝ² = ℝ⁴
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

instance KnownDirectSum ℝ ℝ² where
  {-# INLINE (⊕) #-}
  (⊕) x (ℝ² yz) = ℝ³ $ SIMD.packVector (0,x,y,z) where (y,z) = SIMD.unpackVector yz
  {-# INLINE split #-}
  split (ℝ³ a) = (x, ℝ² $ SIMD.packVector (y,z)) where (_,x,y,z) = SIMD.unpackVector a

instance KnownDirectSum ℝ² ℝ where
  {-# INLINE (⊕) #-}
  (⊕) (ℝ² xy) z = ℝ³ $ SIMD.packVector (0,x,y,z) where (x,y) = SIMD.unpackVector xy
  {-# INLINE split #-}
  split (ℝ³ a) = (ℝ² $ SIMD.packVector (x,y), z) where (_,x,y,z) = SIMD.unpackVector a

instance KnownDirectSum ℝ ℝ³ where
  {-# INLINE (⊕) #-}
  (⊕) x (ℝ³ yzw) = ℝ⁴ $ SIMD.packVector (x,y,z,w) where (_,y,z,w) = SIMD.unpackVector yzw
  {-# INLINE split #-}
  split (ℝ⁴ a) = (x, ℝ³ $ SIMD.packVector (0,y,z,w)) where (x,y,z,w) = SIMD.unpackVector a

instance KnownDirectSum ℝ³ ℝ where
  {-# INLINE (⊕) #-}
  (⊕) (ℝ³ wxy) z = ℝ⁴ $ SIMD.packVector (w,x,y,z) where (_,w,x,y) = SIMD.unpackVector wxy
  {-# INLINE split #-}
  split (ℝ⁴ a) = (ℝ³ $ SIMD.packVector (0,w,x,y), z) where (w,x,y,z) = SIMD.unpackVector a

instance KnownDirectSum ℝ² ℝ² where
  {-# INLINE (⊕) #-}
  (⊕) (ℝ² xy) (ℝ² zw) = ℝ⁴ $ SIMD.packVector (x,y,z,w)
   where (x,y) = SIMD.unpackVector xy; (z,w) = SIMD.unpackVector zw
  {-# INLINE split #-}
  split (ℝ⁴ a) = (ℝ² $ SIMD.packVector (x,y), ℝ² $ SIMD.packVector (z,w))
   where (x,y,z,w) = SIMD.unpackVector a

instance KnownDirectSum (DirectSum u v) w

