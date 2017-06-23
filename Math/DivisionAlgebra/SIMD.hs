-- |
-- Module      : Math.DivisionAlgebra.SIMD
-- Copyright   : (c) Justus Sagemüller 2017
-- License     : GPL v3
-- 
-- Maintainer  : (@) sagemueller $ geo.uni-koeln.de
-- Stability   : experimental
-- Portability : portable
-- 
{-# LANGUAGE TypeFamilies            #-}
{-# LANGUAGE FlexibleInstances       #-}
{-# LANGUAGE FlexibleContexts        #-}
{-# LANGUAGE ScopedTypeVariables     #-}
{-# LANGUAGE UnicodeSyntax           #-}
{-# LANGUAGE DefaultSignatures       #-}
{-# LANGUAGE UndecidableInstances    #-}

module Math.DivisionAlgebra.SIMD where

import Data.AffineSpace
import Data.VectorSpace hiding (magnitudeSq)
import Data.Basis

import qualified Data.Foldable as Foldable

import qualified Data.Primitive.SIMD as SIMD

import qualified Data.Complex as Prelude

import Test.QuickCheck.Arbitrary


data family Complex a

newtype instance Complex Double = ComplexDouble { getComplexFloat :: SIMD.DoubleX2 }

instance Num (Complex Double) where
  {-# INLINE fromInteger #-}
  fromInteger n = ComplexDouble $ SIMD.packVector (fromInteger n, 0)
  {-# INLINE (+) #-}
  ComplexDouble a + ComplexDouble b = ComplexDouble $ a + b
  {-# INLINE (-) #-}
  ComplexDouble a - ComplexDouble b = ComplexDouble $ a - b
  {-# INLINE (*) #-}
  ComplexDouble a * ComplexDouble b = ComplexDouble
      $ a * SIMD.broadcastVector rb + SIMD.packVector (-ia, ra) * SIMD.broadcastVector ib
   where (ra, ia) = SIMD.unpackVector a
         (rb, ib) = SIMD.unpackVector b
  {-# INLINE negate #-}
  negate (ComplexDouble a) = ComplexDouble $ negate a
  {-# INLINE abs #-}
  abs a = modulus a +: 0
  {-# INLINE signum #-}
  signum n@(ComplexDouble a) = ComplexDouble $ a / SIMD.broadcastVector (modulus n)

instance AdditiveGroup (Complex Double) where
  {-# INLINE zeroV #-}
  zeroV = ComplexDouble SIMD.nullVector
  {-# INLINE (^+^) #-}
  (^+^) = (+)
  {-# INLINE negateV #-}
  negateV = negate

-- ^ Note that the scalars are complex, too!
instance VectorSpace (Complex Double) where
  type Scalar (Complex Double) = Complex Double
  {-# INLINE (*^) #-}
  (*^) = (*)

-- ^ @u <.> v  ≡  ̄u * v@ (i.e., sesquilinear and linear in the /second/ argument).
instance InnerSpace (Complex Double) where
  {-# INLINE (<.>) #-}
  ComplexDouble a <.> ComplexDouble b = ComplexDouble
      $ SIMD.packVector (ra,-ia) * SIMD.broadcastVector rb
      + SIMD.packVector (ia, ra) * SIMD.broadcastVector ib
   where (ra, ia) = SIMD.unpackVector a
         (rb, ib) = SIMD.unpackVector b
  
instance Show (Complex Double) where
  showsPrec p (ComplexDouble a) = showParen (p>6) $ shows r . showSign . shows (abs i)
   where (r,i) = SIMD.unpackVector a
         showSign | i<0        = ("-:"++)
                  | otherwise  = ("+:"++)


infixl 6 +:, -:
class (Num c, Num (RealPart c)) => FieldAlgebra c where
  type RealPart c :: *
  type ImagPart c :: *
  (+:) :: RealPart c -> ImagPart c -> c
  (-:) :: RealPart c -> ImagPart c -> c
  realPart :: c -> RealPart c
  imagPart :: c -> ImagPart c
  modulus :: c -> RealPart c
  modulusSq :: c -> RealPart c
  conjugate :: c -> c

instance FieldAlgebra Double where
  type RealPart Double = Double
  type ImagPart Double = ()
  x+:() = x
  x-:() = x
  realPart = id
  imagPart = const ()
  modulus = abs
  modulusSq = (^2)
  conjugate = id
instance FieldAlgebra (Complex Double) where
  type RealPart (Complex Double) = Double
  type ImagPart (Complex Double) = Double
  x+:y = ComplexDouble $ SIMD.packVector (x,y)
  x-:y = ComplexDouble $ SIMD.packVector (x,-y)
  realPart (ComplexDouble a) = r
   where (r, _) = SIMD.unpackVector a
  imagPart (ComplexDouble a) = i
   where (_, i) = SIMD.unpackVector a
  modulus = sqrt . modulusSq
  modulusSq (ComplexDouble a) = SIMD.sumVector $ a*a
  conjugate = \(ComplexDouble a) -> ComplexDouble $ a * cf
   where cf = SIMD.packVector (1,-1)

instance FieldAlgebra Integer where
  type RealPart Integer = Integer
  type ImagPart Integer = ()
  x+:() = x
  x-:() = x
  realPart = id
  imagPart = const ()
  modulus = abs
  modulusSq = (^2)
  conjugate = id
instance FieldAlgebra Int where
  type RealPart Int = Int
  type ImagPart Int = ()
  x+:() = x
  x-:() = x
  realPart = id
  imagPart = const ()
  modulus = abs
  modulusSq = (^2)
  conjugate = id

fromStdComplex :: (FieldAlgebra c, RealPart c ~ s, ImagPart c ~ s)
                    => Prelude.Complex s -> c
fromStdComplex z = Prelude.realPart z +: Prelude.imagPart z

toStdComplex :: (FieldAlgebra c, RealPart c ~ s, ImagPart c ~ s)
                    => c -> Prelude.Complex s
toStdComplex z = realPart z Prelude.:+ imagPart z

instance Arbitrary (Complex Double) where
  arbitrary = fromStdComplex <$> arbitrary
  shrink = map fromStdComplex . shrink . toStdComplex




class (InnerSpace v, FieldAlgebra (Scalar v)) => UnitarySpace v where
  magnitudeSq :: v -> RealPart (Scalar v)
  magnitude :: v -> RealPart (Scalar v)
  default magnitude :: Floating (RealPart (Scalar v)) => v -> RealPart (Scalar v)
  magnitude = sqrt . magnitudeSq


instance UnitarySpace Double where magnitudeSq = (^2)
instance UnitarySpace Integer where
  magnitudeSq = (^2)
  magnitude = abs
instance UnitarySpace Int where
  magnitudeSq = (^2)
  magnitude = abs
instance ( UnitarySpace a, UnitarySpace b
         , Scalar a ~ Scalar b, Floating (RealPart (Scalar a)) )
    => UnitarySpace (a,b) where
  magnitudeSq (x,y) = magnitudeSq x + magnitudeSq y

instance ( FieldAlgebra (Complex s)
         , InnerSpace (Complex s), Scalar (Complex s) ~ Complex s )
                   => UnitarySpace (Complex s) where
  magnitudeSq = modulusSq
  magnitude = modulus





data ComplexT v = ComplexT !v !v

instance (AdditiveGroup v) => AdditiveGroup (ComplexT v) where
  ComplexT ra ia ^+^ ComplexT rb ib = ComplexT (ra^+^rb) (ia^+^ib)
  zeroV = ComplexT zeroV zeroV
  negateV (ComplexT ra ia) = ComplexT (negateV ra) (negateV ia)

instance ( VectorSpace v
         , FieldAlgebra (Complex (Scalar v))
         , RealPart (Complex (Scalar v)) ~ Scalar v
         , ImagPart (Complex (Scalar v)) ~ Scalar v ) => VectorSpace (ComplexT v) where
  type Scalar (ComplexT v) = Complex (Scalar v)
  μ *^ ComplexT ra ia = ComplexT (rμ*^ra ^-^ iμ*^ia)
                                 (iμ*^ra ^+^ rμ*^ia)
   where rμ = realPart μ
         iμ = imagPart μ

instance ( InnerSpace v, AdditiveGroup (Complex (Scalar v))
         , FieldAlgebra (Complex (Scalar v))
         , RealPart (Complex (Scalar v)) ~ Scalar v
         , ImagPart (Complex (Scalar v)) ~ Scalar v ) => InnerSpace (ComplexT v) where
  ComplexT ra ia <.> ComplexT rb ib = (ra<.>rb + ia<.>ib) +: (ra<.>ib - ia<.>rb)

instance ( UnitarySpace v, AdditiveGroup (Complex (Scalar v))
         , Floating (Scalar v), RealPart (Scalar v) ~ Scalar v
         , FieldAlgebra (Complex (Scalar v))
         , RealPart (Complex (Scalar v)) ~ Scalar v
         , ImagPart (Complex (Scalar v)) ~ Scalar v )
         => UnitarySpace (ComplexT v) where
  magnitudeSq (ComplexT ra ia) = magnitudeSq ra + magnitudeSq ia
