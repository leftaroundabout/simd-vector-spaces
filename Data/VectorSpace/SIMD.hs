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
{-# LANGUAGE PatternSynonyms         #-}
{-# LANGUAGE ViewPatterns            #-}

module Data.VectorSpace.SIMD ( ℝ, pattern ℝ
                             , ℝ², pattern ℝ²
                             , ℝ³, pattern ℝ³
                             , ℝ⁴, pattern ℝ⁴
                             , Complex) where

import Math.DivisionAlgebra.SIMD

import Data.AdditiveGroup
import Data.VectorSpace

import qualified Data.Primitive.SIMD as SIMD

import Test.QuickCheck.Arbitrary

type ℝ = Double
newtype ℝ² = R² SIMD.DoubleX2
newtype ℝ³ = R³ SIMD.DoubleX4
newtype ℝ⁴ = R⁴ SIMD.DoubleX4

#define SIMDVSInstances(V,C,F)     \
instance AdditiveGroup (V) where {  \
  C a ^+^ C b = C (a+b);             \
  zeroV = C 0;                        \
  negateV (C a) = C $ negate a };      \
instance VectorSpace V where {          \
  type Scalar (V) = (F);                 \
  μ *^ C a = C                            \
    $ SIMD.broadcastVector μ * a };        \
instance InnerSpace V where {               \
  C a <.> C b = SIMD.sumVector $ a * b };    \
instance UnitarySpace V where {               \
  magnitudeSq (C a) = SIMD.sumVector $ a*a }

SIMDVSInstances(ℝ², R², ℝ)
SIMDVSInstances(ℝ³, R³, ℝ)
SIMDVSInstances(ℝ⁴, R⁴, ℝ)


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
  (⊕) x y = R² $ SIMD.packVector (x,y)
  {-# INLINE split #-}
  split (R² a) = SIMD.unpackVector a

instance KnownDirectSum ℝ ℝ² where
  {-# INLINE (⊕) #-}
  (⊕) x (R² yz) = R³ $ SIMD.packVector (0,x,y,z) where (y,z) = SIMD.unpackVector yz
  {-# INLINE split #-}
  split (R³ a) = (x, R² $ SIMD.packVector (y,z)) where (_,x,y,z) = SIMD.unpackVector a

instance KnownDirectSum ℝ² ℝ where
  {-# INLINE (⊕) #-}
  (⊕) (R² xy) z = R³ $ SIMD.packVector (0,x,y,z) where (x,y) = SIMD.unpackVector xy
  {-# INLINE split #-}
  split (R³ a) = (R² $ SIMD.packVector (x,y), z) where (_,x,y,z) = SIMD.unpackVector a

instance KnownDirectSum ℝ ℝ³ where
  {-# INLINE (⊕) #-}
  (⊕) x (R³ yzw) = R⁴ $ SIMD.packVector (x,y,z,w) where (_,y,z,w) = SIMD.unpackVector yzw
  {-# INLINE split #-}
  split (R⁴ a) = (x, R³ $ SIMD.packVector (0,y,z,w)) where (x,y,z,w) = SIMD.unpackVector a

instance KnownDirectSum ℝ³ ℝ where
  {-# INLINE (⊕) #-}
  (⊕) (R³ wxy) z = R⁴ $ SIMD.packVector (w,x,y,z) where (_,w,x,y) = SIMD.unpackVector wxy
  {-# INLINE split #-}
  split (R⁴ a) = (R³ $ SIMD.packVector (0,w,x,y), z) where (w,x,y,z) = SIMD.unpackVector a

instance KnownDirectSum ℝ² ℝ² where
  {-# INLINE (⊕) #-}
  (⊕) (R² xy) (R² zw) = R⁴ $ SIMD.packVector (x,y,z,w)
   where (x,y) = SIMD.unpackVector xy; (z,w) = SIMD.unpackVector zw
  {-# INLINE split #-}
  split (R⁴ a) = (R² $ SIMD.packVector (x,y), R² $ SIMD.packVector (z,w))
   where (x,y,z,w) = SIMD.unpackVector a

instance KnownDirectSum (DirectSum u v) w



type family V2 a where
  V2 ℝ = ℝ²

-- | Provided for consistency. @ℝ n@ is equivalent to @n@, restricted to type 'Double'.
pattern ℝ :: ℝ -> ℝ
pattern ℝ x = x

pattern ℝ² :: ℝ -> ℝ -> ℝ²
pattern ℝ² x y <- R² (SIMD.unpackVector -> (x,y))
 where ℝ² x y = R² $ SIMD.packVector (x,y)

pattern ℝ³ :: ℝ -> ℝ -> ℝ -> ℝ³
pattern ℝ³ x y z <- R³ (SIMD.unpackVector -> (_,x,y,z))
 where ℝ³ x y z = R³ $ SIMD.packVector (0,x,y,z)

pattern ℝ⁴ :: ℝ -> ℝ -> ℝ -> ℝ -> ℝ⁴
pattern ℝ⁴ x y z w<- R⁴ (SIMD.unpackVector -> (x,y,z,w))
 where ℝ⁴ x y z w = R⁴ $ SIMD.packVector (x,y,z,w)

instance Show ℝ² where
  showsPrec p (ℝ² x y) = showParen (p>9) $ ("ℝ² "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y
instance Show ℝ³ where
  showsPrec p (ℝ³ x y z) = showParen (p>9) $ ("ℝ³ "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y.(' ':).showsPrec 10 z
instance Show ℝ⁴ where
  showsPrec p (ℝ⁴ x y z w) = showParen (p>9) $ ("ℝ⁴ "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y.(' ':).showsPrec 10 z.(' ':).showsPrec 10 w



instance Arbitrary ℝ² where arbitrary = fmap (R² . SIMD.packVector) arbitrary
instance Arbitrary ℝ³ where arbitrary = fmap (R³ . SIMD.packVector) arbitrary
instance Arbitrary ℝ⁴ where arbitrary = fmap (R⁴ . SIMD.packVector) arbitrary
