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
{-# LANGUAGE UndecidableInstances    #-}
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
                             , Complex
                             , ℂ, pattern ℂ
                             , ℂ², pattern ℂ²
                             , ℂ³, pattern ℂ³
                             , ℂ⁴, pattern ℂ⁴
                             ) where

import Math.DivisionAlgebra.SIMD

import Data.AdditiveGroup
import Data.VectorSpace hiding (magnitudeSq)

import qualified Data.Primitive.SIMD as SIMD

import Test.QuickCheck.Arbitrary

type ℝ = Double
newtype ℝ² = R² SIMD.DoubleX2
newtype ℝ³ = R³ SIMD.DoubleX4
newtype ℝ⁴ = R⁴ SIMD.DoubleX4

type ℂ = Complex ℝ
newtype ℂ² = C² SIMD.DoubleX4

type ℂ³ = ComplexT ℝ³
type ℂ⁴ = ComplexT ℝ⁴

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

instance AdditiveGroup ℂ² where
  C² a ^+^ C² b = C² $ a + b
  zeroV = C² 0
  negateV (C² a) = C² $ negate a
instance VectorSpace ℂ² where
  type Scalar ℂ² = ℂ
  ComplexDouble μ *^ C² bc = C²
      $ SIMD.packVector (rμ, iμ, rμ, iμ) * SIMD.packVector (rb, rb, rc, rc)
      + SIMD.packVector (-iμ, rμ, -iμ, rμ) * SIMD.packVector (ib, ib, ic, ic)
   where (rμ, iμ) = SIMD.unpackVector μ
         (rb, ib, rc, ic) = SIMD.unpackVector bc

-- ^ @ℂ² u v <.> ℂ² w x  ≡  ̄u*w + ̄v*x@
instance InnerSpace ℂ² where
  {-# INLINE (<.>) #-}
  C² uv <.> C² wx = ComplexDouble
      $ SIMD.packVector (ru,-iu) * SIMD.broadcastVector rw
      + SIMD.packVector (iu, ru) * SIMD.broadcastVector iw
      + SIMD.packVector (rv,-iv) * SIMD.broadcastVector rx
      + SIMD.packVector (iv, rv) * SIMD.broadcastVector ix
   where (ru,iu,rv,iv) = SIMD.unpackVector uv
         (rw,iw,rx,ix) = SIMD.unpackVector wx

instance UnitarySpace ℂ² where
  magnitudeSq (C² a) = SIMD.sumVector $ a * a





data DirectSum v w = DirectSum !v !w

infixl 6 ⊕
type family (⊕) v w where
  ℝ ⊕ ℝ = ℝ²
  ℝ ⊕ ℝ² = ℝ³
  ℝ² ⊕ ℝ = ℝ³
  ℝ ⊕ ℝ³ = ℝ⁴
  ℝ³ ⊕ ℝ = ℝ⁴
  ℝ² ⊕ ℝ² = ℝ⁴
  ℂ ⊕ ℂ = ℂ²
  ℂ² ⊕ ℂ = ℂ³
  ℂ ⊕ ℂ² = ℂ³
  ℂ² ⊕ ℂ² = ℂ⁴
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

instance KnownDirectSum ℂ ℂ where
  {-# INLINE (⊕) #-}
  (⊕) (ComplexDouble x) (ComplexDouble y)
         = C² $ SIMD.packVector (rx,ix,ry,iy)
   where (rx,ix) = SIMD.unpackVector x
         (ry,iy) = SIMD.unpackVector y
  {-# INLINE split #-}
  split (C² a) = ( ComplexDouble $ SIMD.packVector (rx,ix)
                 , ComplexDouble $ SIMD.packVector (ry,iy) )
   where (rx,ix,ry,iy) = SIMD.unpackVector a

instance KnownDirectSum ℂ ℂ² where
  {-# INLINE (⊕) #-}
  (⊕) (ComplexDouble x) (C² yz)
      = ComplexT (R³ $ SIMD.packVector (0,rx,ry,rz)) (R³ $ SIMD.packVector (0,ix,iy,iz))
   where (rx,ix) = SIMD.unpackVector x
         (ry,iy,rz,iz) = SIMD.unpackVector yz
  {-# INLINE split #-}
  split (ComplexT (R³ ra) (R³ ia)) = ( ComplexDouble $ SIMD.packVector (rx,ix)
                                     , C² $ SIMD.packVector (ry,iy,rz,iz) )
   where (_,rx,ry,rz) = SIMD.unpackVector ra
         (_,ix,iy,iz) = SIMD.unpackVector ia

instance KnownDirectSum ℂ² ℂ where
  {-# INLINE (⊕) #-}
  (⊕) (C² xy) (ComplexDouble z)
      = ComplexT (R³ $ SIMD.packVector (0,rx,ry,rz)) (R³ $ SIMD.packVector (0,ix,iy,iz))
   where (rx,ix,ry,iy) = SIMD.unpackVector xy
         (rz,iz) = SIMD.unpackVector z
  {-# INLINE split #-}
  split (ComplexT (R³ ra) (R³ ia)) = ( C² $ SIMD.packVector (rx,ix,ry,iy)
                                     , ComplexDouble $ SIMD.packVector (rz,iz) )
   where (_,rx,ry,rz) = SIMD.unpackVector ra
         (_,ix,iy,iz) = SIMD.unpackVector ia

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



-- | Provided for consistency. @ℂ n@ is equivalent to @n@, restricted to type 'Complex' 'Double'.
pattern ℂ :: ℂ -> ℂ
pattern ℂ x = x

pattern ℂ² :: ℂ -> ℂ -> ℂ²
pattern ℂ² x y <- (split -> (x,y))
 where ℂ² (ComplexDouble x) (ComplexDouble y) = C² $ SIMD.packVector (rx,ix,ry,iy)
        where (rx,ix) = SIMD.unpackVector x
              (ry,iy) = SIMD.unpackVector y

splitℂ³ :: ℂ³ -> (ℂ,ℂ,ℂ)
splitℂ³ (ComplexT (R³ ra) (R³ ia)) = ( ComplexDouble $ SIMD.packVector (rx,ix)
                                     , ComplexDouble $ SIMD.packVector (ry,iy)
                                     , ComplexDouble $ SIMD.packVector (rz,iz) )
 where (_,rx,ry,rz) = SIMD.unpackVector ra
       (_,ix,iy,iz) = SIMD.unpackVector ia

pattern ℂ³ :: ℂ -> ℂ -> ℂ -> ℂ³
pattern ℂ³ x y z <- (splitℂ³ -> (x,y,z))
 where ℂ³ x y z = ComplexT (ℝ³ (realPart x) (realPart y) (realPart z))
                           (ℝ³ (imagPart x) (imagPart y) (imagPart z))

splitℂ⁴ :: ℂ⁴ -> (ℂ,ℂ,ℂ,ℂ)
splitℂ⁴ (ComplexT (R⁴ ra) (R⁴ ia)) = ( ComplexDouble $ SIMD.packVector (rx,ix)
                                     , ComplexDouble $ SIMD.packVector (ry,iy)
                                     , ComplexDouble $ SIMD.packVector (rz,iz)
                                     , ComplexDouble $ SIMD.packVector (rw,iw) )
 where (rx,ry,rz,rw) = SIMD.unpackVector ra
       (ix,iy,iz,iw) = SIMD.unpackVector ia

pattern ℂ⁴ :: ℂ -> ℂ -> ℂ -> ℂ -> ℂ⁴
pattern ℂ⁴ x y z w <- (splitℂ⁴ -> (x,y,z,w))
 where ℂ⁴ x y z w = ComplexT (ℝ⁴ (realPart x) (realPart y) (realPart z) (realPart w))
                             (ℝ⁴ (imagPart x) (imagPart y) (imagPart z) (imagPart w))


instance Show ℝ² where
  showsPrec p (ℝ² x y) = showParen (p>9) $ ("ℝ² "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y
instance Show ℝ³ where
  showsPrec p (ℝ³ x y z) = showParen (p>9) $ ("ℝ³ "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y.(' ':).showsPrec 10 z
instance Show ℝ⁴ where
  showsPrec p (ℝ⁴ x y z w) = showParen (p>9) $ ("ℝ⁴ "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y.(' ':).showsPrec 10 z.(' ':).showsPrec 10 w
instance Show ℂ² where
  showsPrec p (ℂ² x y) = showParen (p>9) $ ("ℂ² "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y
instance Show ℂ³ where
  showsPrec p (ℂ³ x y z) = showParen (p>9) $ ("ℂ³ "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y.(' ':).showsPrec 10 z
instance Show ℂ⁴ where
  showsPrec p (ℂ⁴ x y z w) = showParen (p>9) $ ("ℂ⁴ "++)
     . showsPrec 10 x.(' ':).showsPrec 10 y.(' ':).showsPrec 10 z.(' ':).showsPrec 10 w



instance Arbitrary ℝ² where arbitrary = fmap (R² . SIMD.packVector) arbitrary
instance Arbitrary ℝ³ where arbitrary = fmap (R³ . SIMD.packVector) arbitrary
instance Arbitrary ℝ⁴ where arbitrary = fmap (R⁴ . SIMD.packVector) arbitrary
instance Arbitrary ℂ² where arbitrary = fmap (C² . SIMD.packVector) arbitrary

instance (Arbitrary v, Arbitrary w) => Arbitrary (DirectSum v w) where
  arbitrary = uncurry DirectSum <$> arbitrary
  shrink (DirectSum v w) = DirectSum <$> shrink v <*> shrink w

instance Arbitrary v => Arbitrary (ComplexT v) where
  arbitrary = uncurry ComplexT <$> arbitrary
  shrink (ComplexT r i) = ComplexT <$> shrink r <*> shrink i
